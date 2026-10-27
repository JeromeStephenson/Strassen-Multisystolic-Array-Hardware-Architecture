//=============================================================================
// smm1_core_dbuf.sv
//
// Double-buffered reusable one-level Strassen datapath core.
//
// This is the double-buffered version of smm1_core.sv.
//
// It implements one level of Strassen recursion:
//
//   1. strassen_pre forms T1..T7 and S1..S7.
//   2. smm1_cores_dbuf computes Q1..Q7 using seven double-buffered baseline
//      systolic cores.
//   3. strassen_post forms C11..C22 from Q1..Q7.
//
// Double-buffer behavior:
//
//   - b_valid / b_row_idx / b_row load the next B/S matrix into the shadow
//     buffers through the Strassen pre-addition network.
//   - swap exchanges active and shadow buffers for all seven lowest cores.
//   - a_valid / a_row stream the current A/T matrix through the active buffers.
//
// This module contains no top-level FSM. It is intended to be instantiated
// inside higher-level tops such as:
//
//   - smm2_top_dbuf.sv
//   - smm2_top_gemm.sv
//
// For SMM2 subcores:
//
//   BASE_N = 8
//   IN_W   = 17
//   MID_W  = 18
//   ACC_W  = 48
//   OUT_W  = 50
//=============================================================================

(* keep_hierarchy = "yes" *)
module smm1_core_dbuf #(
    parameter int BASE_N = 8,
    parameter int IN_W   = 16,
    parameter int MID_W  = IN_W + 1,
    parameter int ACC_W  = 48,
    parameter int OUT_W  = ACC_W + 2
) (
    input  logic                                   clk,
    input  logic                                   rst_n,

    //-------------------------------------------------------------------------
    // B-side shadow-buffer load input
    //
    // b_valid  : valid B row entering this SMM1 core
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
    // core_load_en  : load enable entering the seven lowest double-buffered
    //                 MM cores
    // core_load_row : load row index entering the seven lowest MM cores
    //
    // These are exposed so that a higher-level SMM2 top can count when all
    // B/S rows have finished loading through the recursive pre-addition logic.
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
    // Internal signals
    //-------------------------------------------------------------------------
    logic                                    t_valid;
    logic                                    s_valid;

    logic signed [BASE_N*MID_W-1:0]          t1;
    logic signed [BASE_N*MID_W-1:0]          t2;
    logic signed [BASE_N*MID_W-1:0]          t3;
    logic signed [BASE_N*MID_W-1:0]          t4;
    logic signed [BASE_N*MID_W-1:0]          t5;
    logic signed [BASE_N*MID_W-1:0]          t6;
    logic signed [BASE_N*MID_W-1:0]          t7;

    logic signed [BASE_N*MID_W-1:0]          s1;
    logic signed [BASE_N*MID_W-1:0]          s2;
    logic signed [BASE_N*MID_W-1:0]          s3;
    logic signed [BASE_N*MID_W-1:0]          s4;
    logic signed [BASE_N*MID_W-1:0]          s5;
    logic signed [BASE_N*MID_W-1:0]          s6;
    logic signed [BASE_N*MID_W-1:0]          s7;

    logic                                    q_valid;

    logic signed [BASE_N*ACC_W-1:0]          q1;
    logic signed [BASE_N*ACC_W-1:0]          q2;
    logic signed [BASE_N*ACC_W-1:0]          q3;
    logic signed [BASE_N*ACC_W-1:0]          q4;
    logic signed [BASE_N*ACC_W-1:0]          q5;
    logic signed [BASE_N*ACC_W-1:0]          q6;
    logic signed [BASE_N*ACC_W-1:0]          q7;

    logic                                    c_valid_post;

    logic signed [BASE_N*OUT_W-1:0]          c11_row;
    logic signed [BASE_N*OUT_W-1:0]          c12_row;
    logic signed [BASE_N*OUT_W-1:0]          c21_row;
    logic signed [BASE_N*OUT_W-1:0]          c22_row;

    //-------------------------------------------------------------------------
    // Delayed load row
    //
    // b_valid and b_row_idx arrive on the same cycle.
    // The strassen_pre module adds one pipeline cycle before S rows are valid.
    // Therefore the load row index must also be delayed by one cycle.
    //-------------------------------------------------------------------------
    logic [$clog2(BASE_N)-1:0] load_row_d;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            load_row_d <= '0;
        end else begin
            if (b_valid) begin
                load_row_d <= b_row_idx;
            end
        end
    end

    //-------------------------------------------------------------------------
    // Lowest-level load signals
    //
    // s_valid is b_valid delayed by one cycle through strassen_pre.
    // load_row_d is b_row_idx delayed by one cycle.
    //
    // Important:
    //   The top-level controller must never assert swap while core_load_en
    //   is active.
    //-------------------------------------------------------------------------
    assign core_load_en  = s_valid;
    assign core_load_row = load_row_d;

    //-------------------------------------------------------------------------
    // Strassen pre-addition network
    //
    // For this SMM1 core, G = 1.
    //-------------------------------------------------------------------------
    strassen_pre #(
        .G  (1),
        .N  (BASE_N),
        .W  (IN_W),
        .MW (MID_W)
    ) u_strassen_pre (
        .clk      (clk),
        .rst_n    (rst_n),

        .a_valid  (a_valid),
        .b_valid  (b_valid),

        .a_row    (a_row),
        .b_row    (b_row),

        .t1       (t1),
        .t2       (t2),
        .t3       (t3),
        .t4       (t4),
        .t5       (t5),
        .t6       (t6),
        .t7       (t7),

        .s1       (s1),
        .s2       (s2),
        .s3       (s3),
        .s4       (s4),
        .s5       (s5),
        .s6       (s6),
        .s7       (s7),

        .t_valid  (t_valid),
        .s_valid  (s_valid)
    );

    //-------------------------------------------------------------------------
    // Seven double-buffered Strassen MM cores
    //
    // Q1 = T1 * S1
    // Q2 = T2 * S2
    // Q3 = T3 * S3
    // Q4 = T4 * S4
    // Q5 = T5 * S5
    // Q6 = T6 * S6
    // Q7 = T7 * S7
    //-------------------------------------------------------------------------
    smm1_cores_dbuf #(
        .N     (BASE_N),
        .MW    (MID_W),
        .ACC_W (ACC_W)
    ) u_smm1_cores_dbuf (
        .clk         (clk),
        .rst_n       (rst_n),

        .load_en     (core_load_en),
        .load_row    (core_load_row),

        .s1          (s1),
        .s2          (s2),
        .s3          (s3),
        .s4          (s4),
        .s5          (s5),
        .s6          (s6),
        .s7          (s7),

        .swap        (swap),
        .active_bank (active_bank),

        .t_valid     (t_valid),

        .t1          (t1),
        .t2          (t2),
        .t3          (t3),
        .t4          (t4),
        .t5          (t5),
        .t6          (t6),
        .t7          (t7),

        .q_valid     (q_valid),

        .q1          (q1),
        .q2          (q2),
        .q3          (q3),
        .q4          (q4),
        .q5          (q5),
        .q6          (q6),
        .q7          (q7)
    );

    //-------------------------------------------------------------------------
    // Strassen post-addition network
    //
    // C11 = Q1 + Q4 - Q5 + Q7
    // C12 = Q3 + Q5
    // C21 = Q2 + Q4
    // C22 = Q1 - Q2 + Q3 + Q6
    //-------------------------------------------------------------------------
    strassen_post #(
        .G  (1),
        .N  (BASE_N),
        .QW (ACC_W),
        .CW (OUT_W)
    ) u_strassen_post (
        .clk       (clk),
        .rst_n     (rst_n),

        .in_valid  (q_valid),

        .q1        (q1),
        .q2        (q2),
        .q3        (q3),
        .q4        (q4),
        .q5        (q5),
        .q6        (q6),
        .q7        (q7),

        .c11       (c11_row),
        .c12       (c12_row),
        .c21       (c21_row),
        .c22       (c22_row),

        .out_valid (c_valid_post)
    );

    //-------------------------------------------------------------------------
    // Output C stream
    //
    // Same packing as smm1_core:
    //
    //   c_row = {C22, C21, C12, C11}
    //-------------------------------------------------------------------------
    assign c_valid = c_valid_post;
    assign c_row   = {c22_row, c21_row, c12_row, c11_row};

endmodule
