//=============================================================================
// smm2_top.sv
//
// Top-level SMM2 datapath and controller for one 32x32 matrix tile.
//
// Implements two levels of Strassen recursion:
//
//   Top level:
//     - strassen_pre with G = 4 forms seven T/S block vectors.
//     - Each T/S block vector enters one smm1_core.
//     - Each smm1_core performs the second Strassen recursion level.
//     - strassen_post with G = 4 combines the seven Q block results.
//
// Operation:
//   1. Assert start.
//   2. Feed 8 top-level B row vectors while b_ready is high.
//   3. Feed 8 top-level A row vectors while a_ready is high.
//   4. Receive 8 top-level C row vectors on c_valid / c_row.
//
// Recursive input packing:
//
//   For SMM2, each top-level A/B row vector contains rows from all 16
//   lowest-level subblocks:
//
//     row = {
//       quadrant22 packed row,
//       quadrant21 packed row,
//       quadrant12 packed row,
//       quadrant11 packed row
//     }
//
//   Each quadrant packed row is itself:
//
//     {
//       subquadrant22 row,
//       subquadrant21 row,
//       subquadrant12 row,
//       subquadrant11 row
//     }
//
//   Therefore total input elements per top-level row:
//
//     4 * 4 * 8 = 128 elements
//
//   With DATA_W = 16:
//
//     128 * 16 = 2048 bits
//
// Output packing:
//
//   c_row = {
//     C22 top quadrant packed row,
//     C21 top quadrant packed row,
//     C12 top quadrant packed row,
//     C11 top quadrant packed row
//   }
//=============================================================================

(* keep_hierarchy = "yes" *)
module smm2_top #(
    parameter int TOP_N      = 32,
    parameter int SMM_R      = 2,
    parameter int BLOCK_N    = TOP_N / (1 << SMM_R), // 8 for r=2
    parameter int TOP_G      = 1 << (2*(SMM_R-1)),   // 4 for r=2
    parameter int DATA_W     = 16,

    // Top-level pre-addition bit growth
    parameter int TOP_MID_W  = DATA_W + 1,           // 17

    // SMM1 subcore input width is the top-level Strassen output width
    parameter int SUB_IN_W   = TOP_MID_W,            // 17

    // SMM1 subcore internal pre-addition bit growth
    parameter int SUB_MID_W  = SUB_IN_W + 1,         // 18

    parameter int ACC_W      = 48,

    // Output width from each SMM1 subcore
    parameter int SUB_OUT_W  = ACC_W + 2,            // 50

    // Final SMM2 output width after top-level post-addition
    parameter int OUT_W      = SUB_OUT_W + 2         // 52
) (
    input  logic                                                clk,
    input  logic                                                rst_n,

    input  logic                                                start,
    output logic                                                busy,
    output logic                                                done,

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
    localparam int CNT_W            = $clog2(BLOCK_N) + 1;
    localparam int ROW_IDX_W        = (BLOCK_N > 1) ? $clog2(BLOCK_N) : 1;

    localparam int TOP_QUAD_ELEMS   = TOP_G * BLOCK_N;

    localparam int TOP_IN_BITS      = 4 * TOP_QUAD_ELEMS * DATA_W;
    localparam int TOP_T_BITS       = TOP_QUAD_ELEMS * TOP_MID_W;
    localparam int SUB_IN_BITS      = 4 * BLOCK_N * SUB_IN_W;
    localparam int SUB_OUT_BITS     = 4 * BLOCK_N * SUB_OUT_W;
    localparam int TOP_OUT_BITS     = TOP_QUAD_ELEMS * OUT_W;

    //-------------------------------------------------------------------------
    // FSM
    //-------------------------------------------------------------------------
    typedef enum logic [1:0] {
        ST_IDLE,
        ST_LOAD_B,
        ST_COMPUTE,
        ST_DONE
    } state_e;

    state_e state, state_next;

    //-------------------------------------------------------------------------
    // Counters
    //-------------------------------------------------------------------------
    logic [CNT_W-1:0] b_cnt;
    logic [CNT_W-1:0] load_cnt;
    logic [CNT_W-1:0] a_cnt;
    logic [CNT_W-1:0] c_cnt;

    //-------------------------------------------------------------------------
    // Top-level delayed load row
    //
    // Top pre-addition delays b_valid by one cycle. Therefore the row index
    // must also be delayed by one cycle before entering each smm1_core.
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
    // Arrays for seven SMM1 subcores
    //-------------------------------------------------------------------------
    logic signed [TOP_T_BITS-1:0]  top_t [7];
    logic signed [TOP_T_BITS-1:0]  top_s [7];

    logic [6:0]                    core_load_en;
    logic [ROW_IDX_W-1:0]          core_load_row [7];

    logic [6:0]                    core_c_valid;
    logic signed [SUB_OUT_BITS-1:0] core_c [7];

    logic                          q_valid;

    //-------------------------------------------------------------------------
    // Top-level post-addition signals
    //-------------------------------------------------------------------------
    logic                          c_valid_post;

    logic signed [TOP_OUT_BITS-1:0] c11_top;
    logic signed [TOP_OUT_BITS-1:0] c12_top;
    logic signed [TOP_OUT_BITS-1:0] c21_top;
    logic signed [TOP_OUT_BITS-1:0] c22_top;

    //-------------------------------------------------------------------------
    // Input handshakes
    //-------------------------------------------------------------------------
    assign b_ready = (state == ST_LOAD_B) && (b_cnt < CNT_W'(BLOCK_N));
    assign a_ready = (state == ST_COMPUTE) && (a_cnt < CNT_W'(BLOCK_N));

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
            end
        end
    end

    //-------------------------------------------------------------------------
    // Counters
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            b_cnt    <= '0;
            load_cnt <= '0;
            a_cnt    <= '0;
            c_cnt    <= '0;
        end else begin
            case (state)
                ST_IDLE: begin
                    if (start) begin
                        b_cnt    <= '0;
                        load_cnt <= '0;
                        a_cnt    <= '0;
                        c_cnt    <= '0;
                    end
                end

                ST_LOAD_B: begin
                    if (pre_b_valid) begin
                        b_cnt <= b_cnt + CNT_W'(1);
                    end

                    // Count the actual lowest-level load enables from one
                    // representative SMM1 subcore. All seven receive the same
                    // valid sequence, so core 0 is sufficient for counting.
                    if (core_load_en[0] && (state == ST_LOAD_B)) begin
                        load_cnt <= load_cnt + CNT_W'(1);
                    end
                end

                ST_COMPUTE: begin
                    if (pre_a_valid) begin
                        a_cnt <= a_cnt + CNT_W'(1);
                    end

                    if (c_valid) begin
                        c_cnt <= c_cnt + CNT_W'(1);
                    end
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
                    state_next = ST_LOAD_B;
                end
            end

            ST_LOAD_B: begin
                if (load_cnt == CNT_W'(BLOCK_N)) begin
                    state_next = ST_COMPUTE;
                end
            end

            ST_COMPUTE: begin
                if (c_cnt == CNT_W'(BLOCK_N)) begin
                    state_next = ST_DONE;
                end
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

    assign busy = (state == ST_LOAD_B) || (state == ST_COMPUTE);
    assign done = (state == ST_DONE);

    //-------------------------------------------------------------------------
    // Top-level Strassen pre-addition network
    //
    // For SMM2:
    //
    //   G = 4
    //   N = 8
    //   W = 16
    //   MW = 17
    //
    // Input:
    //
    //   a_row / b_row width = 4 * 4 * 8 * 16 = 2048 bits
    //
    // Output:
    //
    //   each T/S vector width = 4 * 8 * 17 = 544 bits
    //
    // Each T/S vector becomes the input row to one smm1_core.
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
    // Seven SMM1 subcores
    //
    // Each smm1_core performs the second level of Strassen recursion.
    //
    // For each subcore:
    //
    //   BASE_N = 8
    //   IN_W   = 17
    //   MID_W  = 18
    //   ACC_W  = 48
    //   OUT_W  = 50
    //
    // Input row width:
    //
    //   4 * 8 * 17 = 544 bits
    //
    // Output row width:
    //
    //   4 * 8 * 50 = 1600 bits
    //-------------------------------------------------------------------------
    for (genvar gi = 0; gi < 7; gi++) begin : g_smm1_cores
        smm1_core #(
            .BASE_N (BLOCK_N),
            .IN_W   (SUB_IN_W),
            .MID_W  (SUB_MID_W),
            .ACC_W  (ACC_W),
            .OUT_W  (SUB_OUT_W)
        ) u_smm1_core (
            .clk           (clk),
            .rst_n         (rst_n),

            .b_valid       (top_s_valid),
            .b_row_idx     (top_load_row_d),
            .b_row         (top_s[gi]),

            .a_valid       (top_t_valid),
            .a_row         (top_t[gi]),

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
    //
    // For SMM2:
    //
    //   G = 4
    //   N = 8
    //   QW = 50
    //   CW = 52
    //
    // Each Q input is one 1600-bit packed block-product vector from one
    // smm1_core.
    //
    // The post-addition network forms the final top-level C quadrants.
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
    //
    // Same recursive packing:
    //
    //   c_row = {C22, C21, C12, C11}
    //-------------------------------------------------------------------------
    assign c_valid = c_valid_post & (state == ST_COMPUTE);
    assign c_row   = {c22_top, c21_top, c12_top, c11_top};

endmodule
