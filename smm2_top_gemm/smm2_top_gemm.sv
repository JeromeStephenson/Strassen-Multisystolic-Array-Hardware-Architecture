//=============================================================================
// smm2_top_gemm.sv
//
// Paper-style overlapped multi-tile SMM2 top with double-buffered B/S loading
// and current-A / next-B overlap.
//
// For num_tiles = T:
//   B stream: B_tile_0 rows 0..7, B_tile_1 rows 0..7, ... B_tile_T-1 rows 0..7
//   A stream: A_tile_0 rows 0..7, A_tile_1 rows 0..7, ... A_tile_T-1 rows 0..7
//   C stream: C_tile_0 rows 0..7, C_tile_1 rows 0..7, ... C_tile_T-1 rows 0..7
//
// Paper-style behavior:
//   - B_tile_0 is first loaded into the shadow buffer through the recursive
//     Strassen pre-addition network.
//   - Buffers are swapped.
//   - A_tile_0 is multiplied using the active buffer.
//   - While A_tile_0 is being multiplied, B_tile_1 is loaded into the shadow
//     buffer if available.
//   - After C_tile_0 is drained and B_tile_1 is loaded, buffers are swapped.
//   - The process repeats.
//
// Notes:
//   - num_tiles must be stable when start is asserted.
//   - start should be a one-cycle pulse.
//   - c_valid has no backpressure in this bring-up version.
//=============================================================================

(* keep_hierarchy = "yes" *)
module smm2_top_gemm #(
    parameter int TOP_N        = 32,
    parameter int SMM_R        = 2,
    parameter int BLOCK_N      = TOP_N / (1 << SMM_R), // 8 for r=2
    parameter int TOP_G        = 1 << (2*(SMM_R-1)),   // 4 for r=2
    parameter int DATA_W       = 16,

    parameter int TOP_MID_W    = DATA_W + 1,           // 17
    parameter int SUB_IN_W     = TOP_MID_W,            // 17
    parameter int SUB_MID_W    = SUB_IN_W + 1,         // 18

    parameter int ACC_W        = 48,
    parameter int SUB_OUT_W    = ACC_W + 2,            // 50
    parameter int OUT_W        = SUB_OUT_W + 2,        // 52

    parameter int TILE_CNT_W   = 16
) (
    input  logic                                                clk,
    input  logic                                                rst_n,

    input  logic                                                start,
    input  logic [TILE_CNT_W-1:0]                               num_tiles,

    output logic                                                busy,
    output logic                                                done,

    output logic                                                active_bank,

    output logic                                                b_ready,
    input  logic                                                b_valid,
    input  logic signed [4*TOP_G*BLOCK_N*DATA_W-1:0]            b_row,

    output logic                                                a_ready,
    input  logic                                                a_valid,
    input  logic signed [4*TOP_G*BLOCK_N*DATA_W-1:0]            a_row,

    output logic                                                c_valid,
    output logic signed [4*TOP_G*BLOCK_N*OUT_W-1:0]             c_row
);

    //-------------------------------------------------------------------------
    // Local parameters
    //-------------------------------------------------------------------------
    localparam int CNT_W               = $clog2(BLOCK_N) + 1;
    localparam int ROW_IDX_W           = (BLOCK_N > 1) ? $clog2(BLOCK_N) : 1;

    localparam int TOP_QUAD_ELEMS      = TOP_G * BLOCK_N;              // 32
    localparam int TOP_T_BITS          = TOP_QUAD_ELEMS * TOP_MID_W;   // 544
    localparam int CORE_OUT_BITS       = 4 * BLOCK_N * SUB_OUT_W;      // 1600
    localparam int TOP_QUAD_OUT_BITS   = TOP_QUAD_ELEMS * OUT_W;       // 1664

    //-------------------------------------------------------------------------
    // FSM
    //-------------------------------------------------------------------------
    typedef enum logic [2:0] {
        ST_IDLE,
        ST_LOAD_FIRST_B,
        ST_SWAP_FIRST,
        ST_COMPUTE_TILE,
        ST_SWAP_NEXT,
        ST_DONE
    } state_e;

    state_e state, state_next;

    //-------------------------------------------------------------------------
    // Tile and row counters
    //-------------------------------------------------------------------------
    logic [TILE_CNT_W-1:0] tile_idx;
    logic [CNT_W-1:0] b_cnt;
    logic [CNT_W-1:0] s_load_cnt;
    logic [CNT_W-1:0] a_cnt;
    logic [CNT_W-1:0] c_cnt;

    //-------------------------------------------------------------------------
    // Top-level delayed load row
    //
    // Top pre-addition delays b_valid by one cycle. Therefore the row index
    // must also be delayed by one cycle before entering each smm1_core_dbuf.
    //-------------------------------------------------------------------------
    logic [ROW_IDX_W-1:0] top_load_row_d;

    //-------------------------------------------------------------------------
    // Top-level pre-addition signals
    //-------------------------------------------------------------------------
    logic pre_a_valid;
    logic pre_b_valid;

    logic top_t_valid;
    logic top_s_valid;

    logic signed [TOP_T_BITS-1:0] top_t1;
    logic signed [TOP_T_BITS-1:0] top_t2;
    logic signed [TOP_T_BITS-1:0] top_t3;
    logic signed [TOP_T_BITS-1:0] top_t4;
    logic signed [TOP_T_BITS-1:0] top_t5;
    logic signed [TOP_T_BITS-1:0] top_t6;
    logic signed [TOP_T_BITS-1:0] top_t7;

    logic signed [TOP_T_BITS-1:0] top_s1;
    logic signed [TOP_T_BITS-1:0] top_s2;
    logic signed [TOP_T_BITS-1:0] top_s3;
    logic signed [TOP_T_BITS-1:0] top_s4;
    logic signed [TOP_T_BITS-1:0] top_s5;
    logic signed [TOP_T_BITS-1:0] top_s6;
    logic signed [TOP_T_BITS-1:0] top_s7;

    //-------------------------------------------------------------------------
    // Arrays for seven double-buffered SMM1 subcores
    //-------------------------------------------------------------------------
    logic signed [TOP_T_BITS-1:0]   top_t [7];
    logic signed [TOP_T_BITS-1:0]   top_s [7];

    logic [6:0]                     core_load_en;
    logic [ROW_IDX_W-1:0]           core_load_row [7];

    logic [6:0]                     core_active_bank;
    logic [6:0]                     core_c_valid;
    logic signed [CORE_OUT_BITS-1:0] core_c [7];

    logic                           q_valid;

    //-------------------------------------------------------------------------
    // Double-buffer control
    //-------------------------------------------------------------------------
    logic swap;

    assign swap        = (state == ST_SWAP_FIRST) || (state == ST_SWAP_NEXT);
    assign active_bank = core_active_bank[0];

    //-------------------------------------------------------------------------
    // Top-level post-addition signals
    //-------------------------------------------------------------------------
    logic                           c_valid_post;

    logic signed [TOP_QUAD_OUT_BITS-1:0] c11_top;
    logic signed [TOP_QUAD_OUT_BITS-1:0] c12_top;
    logic signed [TOP_QUAD_OUT_BITS-1:0] c21_top;
    logic signed [TOP_QUAD_OUT_BITS-1:0] c22_top;

    //-------------------------------------------------------------------------
    // Tile control helpers
    //-------------------------------------------------------------------------
    logic need_next_b;
    logic last_tile;

    assign need_next_b = ((tile_idx + TILE_CNT_W'(1)) < num_tiles);
    assign last_tile   = (num_tiles == '0)
                       ? 1'b0
                       : (tile_idx == (num_tiles - TILE_CNT_W'(1)));

    //-------------------------------------------------------------------------
    // Input handshakes
    //-------------------------------------------------------------------------
    assign b_ready = ((state == ST_LOAD_FIRST_B) ||
                    ((state == ST_COMPUTE_TILE) && need_next_b))
                && (b_cnt < CNT_W'(BLOCK_N));

    assign a_ready = (state == ST_COMPUTE_TILE)
                && (a_cnt < CNT_W'(BLOCK_N));

    assign pre_b_valid = b_ready & b_valid;
    assign pre_a_valid = a_ready & a_valid;

    //-------------------------------------------------------------------------
    // Top-level load row delay
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            top_load_row_d <= '0;
        end else begin
            if (pre_b_valid) begin
                top_load_row_d <= b_cnt[ROW_IDX_W-1:0];
            end else if ((state == ST_SWAP_FIRST) || (state == ST_SWAP_NEXT)) begin
                top_load_row_d <= '0;
            end
        end
    end

    //-------------------------------------------------------------------------
    // Tile completion condition
    //-------------------------------------------------------------------------
    logic current_tile_c_done;
    logic next_b_loaded_or_not_needed;

    assign current_tile_c_done           = (c_cnt == CNT_W'(BLOCK_N));
    assign next_b_loaded_or_not_needed   = last_tile
                                        || (s_load_cnt == CNT_W'(BLOCK_N));

    //-------------------------------------------------------------------------
    // Counters
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            tile_idx   <= '0;
            b_cnt      <= '0;
            s_load_cnt <= '0;
            a_cnt      <= '0;
            c_cnt      <= '0;
        end else begin
            case (state)
                ST_IDLE: begin
                    if (start) begin
                        tile_idx   <= '0;
                        b_cnt      <= '0;
                        s_load_cnt <= '0;
                        a_cnt      <= '0;
                        c_cnt      <= '0;
                    end
                end

                ST_LOAD_FIRST_B: begin
                    if (pre_b_valid) begin
                        b_cnt <= b_cnt + CNT_W'(1);
                    end

                    // Count the actual lowest-level load enables from one
                    // representative SMM1 subcore. core_load_en[0] already
                    // accounts for the two-level Strassen pre-addition latency.
                    if (core_load_en[0]) begin
                        s_load_cnt <= s_load_cnt + CNT_W'(1);
                    end
                end

                ST_SWAP_FIRST: begin
                    tile_idx   <= '0;
                    b_cnt      <= '0;
                    s_load_cnt <= '0;
                    a_cnt      <= '0;
                    c_cnt      <= '0;
                end

                ST_COMPUTE_TILE: begin
                    if (pre_a_valid) begin
                        a_cnt <= a_cnt + CNT_W'(1);
                    end

                    if (pre_b_valid) begin
                        b_cnt <= b_cnt + CNT_W'(1);
                    end

                    if (core_load_en[0]) begin
                        s_load_cnt <= s_load_cnt + CNT_W'(1);
                    end

                    if (c_valid) begin
                        c_cnt <= c_cnt + CNT_W'(1);
                    end
                end

                ST_SWAP_NEXT: begin
                    tile_idx   <= tile_idx + TILE_CNT_W'(1);
                    b_cnt      <= '0;
                    s_load_cnt <= '0;
                    a_cnt      <= '0;
                    c_cnt      <= '0;
                end

                default: ;
            endcase
        end
    end

    //-------------------------------------------------------------------------
    // FSM next-state logic
    //-------------------------------------------------------------------------
    always_comb begin
        state_next = state;

        case (state)
            ST_IDLE: begin
                if (start) begin
                    if (num_tiles == '0) state_next = ST_DONE;
                    else                 state_next = ST_LOAD_FIRST_B;
                end
            end

            ST_LOAD_FIRST_B: begin
                if (s_load_cnt == CNT_W'(BLOCK_N)) begin
                    state_next = ST_SWAP_FIRST;
                end
            end

            ST_SWAP_FIRST: begin
                state_next = ST_COMPUTE_TILE;
            end

            ST_COMPUTE_TILE: begin
                if (current_tile_c_done && next_b_loaded_or_not_needed) begin
                    if (last_tile) state_next = ST_DONE;
                    else           state_next = ST_SWAP_NEXT;
                end
            end

            ST_SWAP_NEXT: begin
                state_next = ST_COMPUTE_TILE;
            end

            ST_DONE: begin
                state_next = ST_IDLE;
            end

            default: begin
                state_next = ST_IDLE;
            end
        endcase
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            state <= ST_IDLE;
        end else begin
            state <= state_next;
        end
    end

    assign busy = (state != ST_IDLE) && (state != ST_DONE);
    assign done = (state == ST_DONE);

    //-------------------------------------------------------------------------
    // Top-level Strassen pre-addition network
    //-------------------------------------------------------------------------
    strassen_pre #(
        .G  (TOP_G),
        .N  (BLOCK_N),
        .W  (DATA_W),
        .MW (TOP_MID_W)
    ) u_top_strassen_pre (
        .clk      (clk),
        .rst_n    (rst_n),

        .a_valid  (pre_a_valid),
        .b_valid  (pre_b_valid),

        .a_row    (a_row),
        .b_row    (b_row),

        .t1       (top_t1),
        .t2       (top_t2),
        .t3       (top_t3),
        .t4       (top_t4),
        .t5       (top_t5),
        .t6       (top_t6),
        .t7       (top_t7),

        .s1       (top_s1),
        .s2       (top_s2),
        .s3       (top_s3),
        .s4       (top_s4),
        .s5       (top_s5),
        .s6       (top_s6),
        .s7       (top_s7),

        .t_valid  (top_t_valid),
        .s_valid  (top_s_valid)
    );

    //-------------------------------------------------------------------------
    // Assign top-level T/S vectors into arrays
    //-------------------------------------------------------------------------
    assign top_t[0] = top_t1;
    assign top_t[1] = top_t2;
    assign top_t[2] = top_t3;
    assign top_t[3] = top_t4;
    assign top_t[4] = top_t5;
    assign top_t[5] = top_t6;
    assign top_t[6] = top_t7;

    assign top_s[0] = top_s1;
    assign top_s[1] = top_s2;
    assign top_s[2] = top_s3;
    assign top_s[3] = top_s4;
    assign top_s[4] = top_s5;
    assign top_s[5] = top_s6;
    assign top_s[6] = top_s7;

    //-------------------------------------------------------------------------
    // Seven double-buffered SMM1 subcores
    //-------------------------------------------------------------------------
    for (genvar gi = 0; gi < 7; gi++) begin : g_smm1_cores_dbuf
        smm1_core_dbuf #(
            .BASE_N (BLOCK_N),
            .IN_W   (SUB_IN_W),
            .MID_W  (SUB_MID_W),
            .ACC_W  (ACC_W),
            .OUT_W  (SUB_OUT_W)
        ) u_smm1_core_dbuf (
            .clk           (clk),
            .rst_n         (rst_n),

            .b_valid       (top_s_valid),
            .b_row_idx     (top_load_row_d),
            .b_row         (top_s[gi]),

            .a_valid       (top_t_valid),
            .a_row         (top_t[gi]),

            .swap          (swap),
            .active_bank   (core_active_bank[gi]),

            .core_load_en  (core_load_en[gi]),
            .core_load_row (core_load_row[gi]),

            .c_valid       (core_c_valid[gi]),
            .c_row         (core_c[gi])
        );
    end

    //-------------------------------------------------------------------------
    // All seven SMM1 subcores have identical latency.
    //-------------------------------------------------------------------------
    assign q_valid = &core_c_valid;

    //-------------------------------------------------------------------------
    // Top-level Strassen post-addition network
    //-------------------------------------------------------------------------
    strassen_post #(
        .G  (TOP_G),
        .N  (BLOCK_N),
        .QW (SUB_OUT_W),
        .CW (OUT_W)
    ) u_top_strassen_post (
        .clk       (clk),
        .rst_n     (rst_n),

        .in_valid  (q_valid),

        .q1        (core_c[0]),
        .q2        (core_c[1]),
        .q3        (core_c[2]),
        .q4        (core_c[3]),
        .q5        (core_c[4]),
        .q6        (core_c[5]),
        .q7        (core_c[6]),

        .c11       (c11_top),
        .c12       (c12_top),
        .c21       (c21_top),
        .c22       (c22_top),

        .out_valid (c_valid_post)
    );

    //-------------------------------------------------------------------------
    // Output C stream
    //-------------------------------------------------------------------------
    assign c_valid = c_valid_post & (state == ST_COMPUTE_TILE);
    assign c_row   = {c22_top, c21_top, c12_top, c11_top};

endmodule
