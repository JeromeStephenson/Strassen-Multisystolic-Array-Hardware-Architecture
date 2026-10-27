//=============================================================================
// mm1_core_dbuf.sv
//
// Double-buffered reusable conventional blocked matrix multiplication core.
//
// This is the baseline conventional double-buffered equivalent of
// smm1_core_dbuf.sv.
//
// It multiplies two 16x16 blocks using eight double-buffered 8x8 baseline
// systolic cores:
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
// Then conventional post-addition forms:
//
//   C11 = Q1 + Q3
//   C12 = Q2 + Q4
//   C21 = Q5 + Q7
//   C22 = Q6 + Q8
//
// Double-buffer behavior:
//
//   - b_valid / b_row_idx / b_row load the next B matrix into the shadow
//     buffers.
//   - swap exchanges active and shadow buffers for all eight lowest cores.
//   - a_valid / a_row stream the current A matrix through the active buffers.
//
// This module contains no top-level FSM. It is intended to be instantiated
// inside higher-level tops such as:
//
//   - mm2_top_dbuf.sv
//   - mm2_top_gemm.sv
//
// For MM2 subcores:
//
//   BASE_N = 8
//   IN_W   = 16
//   ACC_W  = 48
//   OUT_W  = 50
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm1_core_dbuf #(
    parameter int BASE_N = 8,
    parameter int IN_W   = 16,
    parameter int ACC_W  = 48,
    parameter int OUT_W  = ACC_W + 2
) (
    input  logic                                   clk,
    input  logic                                   rst_n,

    //-------------------------------------------------------------------------
    // B-side shadow-buffer load input
    //
    // b_valid  : valid B row entering this MM1 core
    // b_row_idx: row index inside the BASE_N x BASE_N lowest blocks
    // b_row    : packed B row containing four subblock quadrants
    //-------------------------------------------------------------------------
    input  logic                                   b_valid,
    input  logic [$clog2(BASE_N)-1:0]              b_row_idx,
    input  logic signed [4*BASE_N*IN_W-1:0]        b_row,

    //-------------------------------------------------------------------------
    // A-side streaming input
    //-------------------------------------------------------------------------
    input  logic                                   a_valid,
    input  logic signed [4*BASE_N*IN_W-1:0]        a_row,

    //-------------------------------------------------------------------------
    // Double-buffer control
    //
    // swap        : one-cycle pulse, exchanges active and shadow buffers
    // active_bank : 0 = bank 0 active, 1 = bank 1 active
    //-------------------------------------------------------------------------
    input  logic                                   swap,
    output logic                                   active_bank,

    //-------------------------------------------------------------------------
    // Internal lowest-level load observation
    //
    // core_load_en  : load enable entering the eight lowest double-buffered
    //                 MM cores
    // core_load_row : load row index entering the eight lowest MM cores
    //
    // For the baseline MM1 core there is no Strassen pre-addition latency,
    // so core_load_en is equal to b_valid.
    //-------------------------------------------------------------------------
    output logic                                   core_load_en,
    output logic [$clog2(BASE_N)-1:0]              core_load_row,

    //-------------------------------------------------------------------------
    // C-side output
    //-------------------------------------------------------------------------
    output logic                                   c_valid,
    output logic signed [4*BASE_N*OUT_W-1:0]       c_row
);

    //-------------------------------------------------------------------------
    // Local parameters
    //-------------------------------------------------------------------------
    localparam int QW = BASE_N * IN_W;
    localparam int CW = OUT_W;

    //-------------------------------------------------------------------------
    // Unpack A quadrants
    //
    // Packing order:
    //
    //   a_row = {A22, A21, A12, A11}
    //
    // Therefore:
    //
    //   A11 is at offset 0
    //   A12 is at offset 1
    //   A21 is at offset 2
    //   A22 is at offset 3
    //-------------------------------------------------------------------------
    logic signed [QW-1:0] a11;
    logic signed [QW-1:0] a12;
    logic signed [QW-1:0] a21;
    logic signed [QW-1:0] a22;

    assign a11 = a_row[0*QW +: QW];
    assign a12 = a_row[1*QW +: QW];
    assign a21 = a_row[2*QW +: QW];
    assign a22 = a_row[3*QW +: QW];

    //-------------------------------------------------------------------------
    // Unpack B quadrants
    //-------------------------------------------------------------------------
    logic signed [QW-1:0] b11;
    logic signed [QW-1:0] b12;
    logic signed [QW-1:0] b21;
    logic signed [QW-1:0] b22;

    assign b11 = b_row[0*QW +: QW];
    assign b12 = b_row[1*QW +: QW];
    assign b21 = b_row[2*QW +: QW];
    assign b22 = b_row[3*QW +: QW];

    //-------------------------------------------------------------------------
    // Baseline load signals
    //
    // No Strassen pre-addition latency exists in the baseline MM1 core.
    //-------------------------------------------------------------------------
    assign core_load_en  = b_valid;
    assign core_load_row = b_row_idx;

    //-------------------------------------------------------------------------
    // Internal valid signals
    //-------------------------------------------------------------------------
    logic q1_valid;
    logic q2_valid;
    logic q3_valid;
    logic q4_valid;
    logic q5_valid;
    logic q6_valid;
    logic q7_valid;
    logic q8_valid;

    logic q_valid;

    //-------------------------------------------------------------------------
    // Active bank status
    //
    // All eight cores receive the same swap signal, so all active_bank
    // outputs should always match. We expose core 0 as the representative
    // status bit.
    //-------------------------------------------------------------------------
    logic [7:0] core_active_bank;

    assign active_bank = core_active_bank[0];

    //-------------------------------------------------------------------------
    // Q outputs from the eight lowest-level cores
    //-------------------------------------------------------------------------
    logic signed [BASE_N*ACC_W-1:0] q1;
    logic signed [BASE_N*ACC_W-1:0] q2;
    logic signed [BASE_N*ACC_W-1:0] q3;
    logic signed [BASE_N*ACC_W-1:0] q4;
    logic signed [BASE_N*ACC_W-1:0] q5;
    logic signed [BASE_N*ACC_W-1:0] q6;
    logic signed [BASE_N*ACC_W-1:0] q7;
    logic signed [BASE_N*ACC_W-1:0] q8;

    //-------------------------------------------------------------------------
    // Q1 = A11 * B11
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q1 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b11),

        .swap        (swap),
        .active_bank (core_active_bank[0]),

        .a_valid     (a_valid),
        .a_row       (a11),

        .c_valid     (q1_valid),
        .c_row       (q1)
    );

    //-------------------------------------------------------------------------
    // Q2 = A11 * B12
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q2 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b12),

        .swap        (swap),
        .active_bank (core_active_bank[1]),

        .a_valid     (a_valid),
        .a_row       (a11),

        .c_valid     (q2_valid),
        .c_row       (q2)
    );

    //-------------------------------------------------------------------------
    // Q3 = A12 * B21
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q3 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b21),

        .swap        (swap),
        .active_bank (core_active_bank[2]),

        .a_valid     (a_valid),
        .a_row       (a12),

        .c_valid     (q3_valid),
        .c_row       (q3)
    );

    //-------------------------------------------------------------------------
    // Q4 = A12 * B22
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q4 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b22),

        .swap        (swap),
        .active_bank (core_active_bank[3]),

        .a_valid     (a_valid),
        .a_row       (a12),

        .c_valid     (q4_valid),
        .c_row       (q4)
    );

    //-------------------------------------------------------------------------
    // Q5 = A21 * B11
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q5 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b11),

        .swap        (swap),
        .active_bank (core_active_bank[4]),

        .a_valid     (a_valid),
        .a_row       (a21),

        .c_valid     (q5_valid),
        .c_row       (q5)
    );

    //-------------------------------------------------------------------------
    // Q6 = A21 * B12
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q6 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b12),

        .swap        (swap),
        .active_bank (core_active_bank[5]),

        .a_valid     (a_valid),
        .a_row       (a21),

        .c_valid     (q6_valid),
        .c_row       (q6)
    );

    //-------------------------------------------------------------------------
    // Q7 = A22 * B21
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q7 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b21),

        .swap        (swap),
        .active_bank (core_active_bank[6]),

        .a_valid     (a_valid),
        .a_row       (a22),

        .c_valid     (q7_valid),
        .c_row       (q7)
    );

    //-------------------------------------------------------------------------
    // Q8 = A22 * B22
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X     (BASE_N),
        .Y     (BASE_N),
        .A_W   (IN_W),
        .B_W   (IN_W),
        .ACC_W (ACC_W)
    ) u_q8 (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),
        .b_row       (b22),

        .swap        (swap),
        .active_bank (core_active_bank[7]),

        .a_valid     (a_valid),
        .a_row       (a22),

        .c_valid     (q8_valid),
        .c_row       (q8)
    );

    //-------------------------------------------------------------------------
    // All eight cores have identical latency.
    //-------------------------------------------------------------------------
    assign q_valid = q1_valid & q2_valid & q3_valid & q4_valid
                   & q5_valid & q6_valid & q7_valid & q8_valid;

    //-------------------------------------------------------------------------
    // Conventional post-addition
    //
    //   C11 = Q1 + Q3
    //   C12 = Q2 + Q4
    //   C21 = Q5 + Q7
    //   C22 = Q6 + Q8
    //
    // Two pipeline stages are used to mirror the timing style of the Strassen
    // post-addition network.
    //-------------------------------------------------------------------------
    function automatic logic signed [CW-1:0] sextq(
        input logic signed [ACC_W-1:0] x
    );
        if (CW == ACC_W) begin
            sextq = x;
        end else begin
            sextq = $signed({{(CW-ACC_W){x[ACC_W-1]}}, x});
        end
    endfunction

    logic v1;
    logic c_valid_post;

    logic signed [BASE_N*ACC_W-1:0] q1r;
    logic signed [BASE_N*ACC_W-1:0] q2r;
    logic signed [BASE_N*ACC_W-1:0] q3r;
    logic signed [BASE_N*ACC_W-1:0] q4r;
    logic signed [BASE_N*ACC_W-1:0] q5r;
    logic signed [BASE_N*ACC_W-1:0] q6r;
    logic signed [BASE_N*ACC_W-1:0] q7r;
    logic signed [BASE_N*ACC_W-1:0] q8r;

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

            q1r <= q1;
            q2r <= q2;
            q3r <= q3;
            q4r <= q4;
            q5r <= q5;
            q6r <= q6;
            q7r <= q7;
            q8r <= q8;
        end
    end

    logic signed [BASE_N*CW-1:0] c11_row;
    logic signed [BASE_N*CW-1:0] c12_row;
    logic signed [BASE_N*CW-1:0] c21_row;
    logic signed [BASE_N*CW-1:0] c22_row;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c11_row <= '0;
            c12_row <= '0;
            c21_row <= '0;
            c22_row <= '0;
        end else begin
            for (int i = 0; i < BASE_N; i++) begin
                c11_row[i*CW +: CW] <=
                    sextq(q1r[i*ACC_W +: ACC_W]) +
                    sextq(q3r[i*ACC_W +: ACC_W]);

                c12_row[i*CW +: CW] <=
                    sextq(q2r[i*ACC_W +: ACC_W]) +
                    sextq(q4r[i*ACC_W +: ACC_W]);

                c21_row[i*CW +: CW] <=
                    sextq(q5r[i*ACC_W +: ACC_W]) +
                    sextq(q7r[i*ACC_W +: ACC_W]);

                c22_row[i*CW +: CW] <=
                    sextq(q6r[i*ACC_W +: ACC_W]) +
                    sextq(q8r[i*ACC_W +: ACC_W]);
            end
        end
    end

    //-------------------------------------------------------------------------
    // Output C stream
    //
    // Same packing as mm1_core:
    //
    //   c_row = {C22, C21, C12, C11}
    //-------------------------------------------------------------------------
    assign c_valid = c_valid_post;
    assign c_row   = {c22_row, c21_row, c12_row, c11_row};

endmodule
