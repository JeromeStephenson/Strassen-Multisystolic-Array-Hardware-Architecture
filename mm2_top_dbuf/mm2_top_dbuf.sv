//=============================================================================
// mm2_top_dbuf.sv
//
// Double-buffered baseline conventional two-level multisystolic array top
// (MM2) for one 32x32 tile.
//
// This is the baseline equivalent of smm2_top_dbuf.sv.
//
// Top-level conventional blocked matrix multiplication:
//
//   Q1 = A11 * B11
//   Q2 = A11 * B12
//   Q3 = A12 * B21
//   Q4 = A12 * B22
//   Q5 = A21 * B11
//   Q6 = A21 * B12
//   Q7 = A22 * B21
//   Q8 = A22 * B22
//
//   C11 = Q1 + Q3
//   C12 = Q2 + Q4
//   C21 = Q5 + Q7
//   C22 = Q6 + Q8
//
// Each Qi is computed by one mm1_core_dbuf.
// Each mm1_core_dbuf contains eight double-buffered 8x8 baseline systolic
// cores.
//
// Double-buffer operation for one tile:
//   1. Assert start.
//   2. Feed BLOCK_N top-level B rows while b_ready is high.
//   3. B rows are loaded into the shadow buffers.
//   4. The FSM asserts swap for one cycle.
//   5. Feed BLOCK_N top-level A rows while a_ready is high.
//   6. Receive BLOCK_N top-level C rows on c_valid / c_row.
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm2_top_dbuf #(
    parameter int TOP_N        = 32,
    parameter int SMM_R        = 2,
    parameter int BLOCK_N      = TOP_N / (1 << SMM_R), // 8 for r=2
    parameter int TOP_G        = 1 << (2*(SMM_R-1)),   // 4 for r=2
    parameter int DATA_W       = 16,

    parameter int ACC_W        = 48,

    // Output width from each mm1_core_dbuf
    parameter int CORE_OUT_W   = ACC_W + 2,            // 50

    // Final MM2 output width after top-level post-addition
    parameter int OUT_W        = CORE_OUT_W + 2        // 52
) (
    input  logic                                                clk,
    input  logic                                                rst_n,

    input  logic                                                start,
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
    localparam int TOP_QUAD_IN_BITS    = TOP_QUAD_ELEMS * DATA_W;      // 512
    localparam int CORE_OUT_BITS       = 4 * BLOCK_N * CORE_OUT_W;     // 1600
    localparam int TOP_QUAD_OUT_BITS   = TOP_QUAD_ELEMS * OUT_W;       // 1664

    //-------------------------------------------------------------------------
    // FSM
    //-------------------------------------------------------------------------
    typedef enum logic [2:0] {
        ST_IDLE,
        ST_LOAD_B,
        ST_SWAP,
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
    // Handshakes
    //-------------------------------------------------------------------------
    logic pre_a_valid;
    logic pre_b_valid;

    assign b_ready = (state == ST_LOAD_B) && (b_cnt < CNT_W'(BLOCK_N));
    assign a_ready = (state == ST_COMPUTE) && (a_cnt < CNT_W'(BLOCK_N));

    assign pre_b_valid = b_ready & b_valid;
    assign pre_a_valid = a_ready & a_valid;

    //-------------------------------------------------------------------------
    // Double-buffer control
    //-------------------------------------------------------------------------
    logic swap;

    assign swap        = (state == ST_SWAP);

    //-------------------------------------------------------------------------
    // Unpack top-level A quadrants
    //
    // Packing order:
    //
    //   a_row = {A22_top, A21_top, A12_top, A11_top}
    //
    // Therefore:
    //
    //   A11_top is at offset 0
    //   A12_top is at offset 1
    //   A21_top is at offset 2
    //   A22_top is at offset 3
    //-------------------------------------------------------------------------
    logic signed [TOP_QUAD_IN_BITS-1:0] a11_top;
    logic signed [TOP_QUAD_IN_BITS-1:0] a12_top;
    logic signed [TOP_QUAD_IN_BITS-1:0] a21_top;
    logic signed [TOP_QUAD_IN_BITS-1:0] a22_top;

    assign a11_top = a_row[0*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];
    assign a12_top = a_row[1*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];
    assign a21_top = a_row[2*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];
    assign a22_top = a_row[3*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];

    //-------------------------------------------------------------------------
    // Unpack top-level B quadrants
    //-------------------------------------------------------------------------
    logic signed [TOP_QUAD_IN_BITS-1:0] b11_top;
    logic signed [TOP_QUAD_IN_BITS-1:0] b12_top;
    logic signed [TOP_QUAD_IN_BITS-1:0] b21_top;
    logic signed [TOP_QUAD_IN_BITS-1:0] b22_top;

    assign b11_top = b_row[0*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];
    assign b12_top = b_row[1*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];
    assign b21_top = b_row[2*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];
    assign b22_top = b_row[3*TOP_QUAD_IN_BITS +: TOP_QUAD_IN_BITS];

    //-------------------------------------------------------------------------
    // Eight double-buffered MM1 subcores
    //-------------------------------------------------------------------------
    logic [7:0]                      core_load_en;
    logic [ROW_IDX_W-1:0]            core_load_row [8];

    logic [7:0]                      core_active_bank;
    logic [7:0]                      core_c_valid;
    logic signed [CORE_OUT_BITS-1:0] core_c [8];

    assign active_bank = core_active_bank[0];

    //-------------------------------------------------------------------------
    // Q1 = A11 * B11
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q1 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b11_top),

        .a_valid       (pre_a_valid),
        .a_row         (a11_top),

        .swap          (swap),
        .active_bank   (core_active_bank[0]),

        .core_load_en  (core_load_en[0]),
        .core_load_row (core_load_row[0]),

        .c_valid       (core_c_valid[0]),
        .c_row         (core_c[0])
    );

    //-------------------------------------------------------------------------
    // Q2 = A11 * B12
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q2 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b12_top),

        .a_valid       (pre_a_valid),
        .a_row         (a11_top),

        .swap          (swap),
        .active_bank   (core_active_bank[1]),

        .core_load_en  (core_load_en[1]),
        .core_load_row (core_load_row[1]),

        .c_valid       (core_c_valid[1]),
        .c_row         (core_c[1])
    );

    //-------------------------------------------------------------------------
    // Q3 = A12 * B21
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q3 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b21_top),

        .a_valid       (pre_a_valid),
        .a_row         (a12_top),

        .swap          (swap),
        .active_bank   (core_active_bank[2]),

        .core_load_en  (core_load_en[2]),
        .core_load_row (core_load_row[2]),

        .c_valid       (core_c_valid[2]),
        .c_row         (core_c[2])
    );

    //-------------------------------------------------------------------------
    // Q4 = A12 * B22
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q4 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b22_top),

        .a_valid       (pre_a_valid),
        .a_row         (a12_top),

        .swap          (swap),
        .active_bank   (core_active_bank[3]),

        .core_load_en  (core_load_en[3]),
        .core_load_row (core_load_row[3]),

        .c_valid       (core_c_valid[3]),
        .c_row         (core_c[3])
    );

    //-------------------------------------------------------------------------
    // Q5 = A21 * B11
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q5 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b11_top),

        .a_valid       (pre_a_valid),
        .a_row         (a21_top),

        .swap          (swap),
        .active_bank   (core_active_bank[4]),

        .core_load_en  (core_load_en[4]),
        .core_load_row (core_load_row[4]),

        .c_valid       (core_c_valid[4]),
        .c_row         (core_c[4])
    );

    //-------------------------------------------------------------------------
    // Q6 = A21 * B12
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q6 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b12_top),

        .a_valid       (pre_a_valid),
        .a_row         (a21_top),

        .swap          (swap),
        .active_bank   (core_active_bank[5]),

        .core_load_en  (core_load_en[5]),
        .core_load_row (core_load_row[5]),

        .c_valid       (core_c_valid[5]),
        .c_row         (core_c[5])
    );

    //-------------------------------------------------------------------------
    // Q7 = A22 * B21
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q7 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b21_top),

        .a_valid       (pre_a_valid),
        .a_row         (a22_top),

        .swap          (swap),
        .active_bank   (core_active_bank[6]),

        .core_load_en  (core_load_en[6]),
        .core_load_row (core_load_row[6]),

        .c_valid       (core_c_valid[6]),
        .c_row         (core_c[6])
    );

    //-------------------------------------------------------------------------
    // Q8 = A22 * B22
    //-------------------------------------------------------------------------
    mm1_core_dbuf #(
        .BASE_N (BLOCK_N),
        .IN_W   (DATA_W),
        .ACC_W  (ACC_W),
        .OUT_W  (CORE_OUT_W)
    ) u_q8 (
        .clk           (clk),
        .rst_n         (rst_n),

        .b_valid       (pre_b_valid),
        .b_row_idx     (b_cnt[ROW_IDX_W-1:0]),
        .b_row         (b22_top),

        .a_valid       (pre_a_valid),
        .a_row         (a22_top),

        .swap          (swap),
        .active_bank   (core_active_bank[7]),

        .core_load_en  (core_load_en[7]),
        .core_load_row (core_load_row[7]),

        .c_valid       (core_c_valid[7]),
        .c_row         (core_c[7])
    );

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

                    // Baseline MM1 cores have no Strassen pre-addition latency.
                    // core_load_en is equal to pre_b_valid inside mm1_core_dbuf.
                    if (core_load_en[0]) begin
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
                    state_next = ST_SWAP;
                end
            end

            ST_SWAP: begin
                state_next = ST_COMPUTE;
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

    assign busy = (state == ST_LOAD_B) || (state == ST_SWAP) || (state == ST_COMPUTE);
    assign done = (state == ST_DONE);

    //-------------------------------------------------------------------------
    // All eight MM1 subcores have identical latency.
    //-------------------------------------------------------------------------
    logic q_valid;

    assign q_valid = &core_c_valid;

    //-------------------------------------------------------------------------
    // Top-level conventional post-addition
    //
    //   C11_top = Q1 + Q3
    //   C12_top = Q2 + Q4
    //   C21_top = Q5 + Q7
    //   C22_top = Q6 + Q8
    //
    // Two pipeline stages are used to mirror the timing style of the
    // Strassen post-addition network.
    //-------------------------------------------------------------------------
    function automatic logic signed [OUT_W-1:0] sext_core(
        input logic signed [CORE_OUT_W-1:0] x
    );
        if (OUT_W == CORE_OUT_W) begin
            sext_core = x;
        end else begin
            sext_core = $signed({{(OUT_W-CORE_OUT_W){x[CORE_OUT_W-1]}}, x});
        end
    endfunction

    logic v1;
    logic c_valid_post;

    logic signed [CORE_OUT_BITS-1:0] q1r;
    logic signed [CORE_OUT_BITS-1:0] q2r;
    logic signed [CORE_OUT_BITS-1:0] q3r;
    logic signed [CORE_OUT_BITS-1:0] q4r;
    logic signed [CORE_OUT_BITS-1:0] q5r;
    logic signed [CORE_OUT_BITS-1:0] q6r;
    logic signed [CORE_OUT_BITS-1:0] q7r;
    logic signed [CORE_OUT_BITS-1:0] q8r;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            v1           <= 1'b0;
            c_valid_post <= 1'b0;

            q1r <= '0;
            q2r <= '0;
            q3r <= '0;
            q4r <= '0;
            q5r <= '0;
            q6r <= '0;
            q7r <= '0;
            q8r <= '0;
        end else begin
            v1           <= q_valid;
            c_valid_post <= v1;

            q1r <= core_c[0];
            q2r <= core_c[1];
            q3r <= core_c[2];
            q4r <= core_c[3];
            q5r <= core_c[4];
            q6r <= core_c[5];
            q7r <= core_c[6];
            q8r <= core_c[7];
        end
    end

    logic signed [TOP_QUAD_OUT_BITS-1:0] c11_top;
    logic signed [TOP_QUAD_OUT_BITS-1:0] c12_top;
    logic signed [TOP_QUAD_OUT_BITS-1:0] c21_top;
    logic signed [TOP_QUAD_OUT_BITS-1:0] c22_top;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c11_top <= '0;
            c12_top <= '0;
            c21_top <= '0;
            c22_top <= '0;
        end else begin
            for (int i = 0; i < TOP_QUAD_ELEMS; i++) begin
                c11_top[i*OUT_W +: OUT_W] <=
                    sext_core(q1r[i*CORE_OUT_W +: CORE_OUT_W]) +
                    sext_core(q3r[i*CORE_OUT_W +: CORE_OUT_W]);

                c12_top[i*OUT_W +: OUT_W] <=
                    sext_core(q2r[i*CORE_OUT_W +: CORE_OUT_W]) +
                    sext_core(q4r[i*CORE_OUT_W +: CORE_OUT_W]);

                c21_top[i*OUT_W +: OUT_W] <=
                    sext_core(q5r[i*CORE_OUT_W +: CORE_OUT_W]) +
                    sext_core(q7r[i*CORE_OUT_W +: CORE_OUT_W]);

                c22_top[i*OUT_W +: OUT_W] <=
                    sext_core(q6r[i*CORE_OUT_W +: CORE_OUT_W]) +
                    sext_core(q8r[i*CORE_OUT_W +: CORE_OUT_W]);
            end
        end
    end

    //-------------------------------------------------------------------------
    // Output C stream
    //
    // Same recursive packing as smm2_top_dbuf:
    //
    //   c_row = {C22_top, C21_top, C12_top, C11_top}
    //-------------------------------------------------------------------------
    assign c_valid = c_valid_post & (state == ST_COMPUTE);
    assign c_row   = {c22_top, c21_top, c12_top, c11_top};

endmodule

