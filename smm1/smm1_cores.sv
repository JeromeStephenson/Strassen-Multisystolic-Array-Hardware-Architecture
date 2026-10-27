//=============================================================================
// smm1_cores.sv
//
// Seven-core Strassen MM array for SMM1.
//
// Instantiates seven baseline systolic MM cores to compute:
//
//   Q1 = T1 * S1
//   Q2 = T2 * S2
//   Q3 = T3 * S3
//   Q4 = T4 * S4
//   Q5 = T5 * S5
//   Q6 = T6 * S6
//   Q7 = T7 * S7
//
// Target configuration:
//   TOP_N   = 16
//   BLOCK_N = 8
//   MID_W   = 17
//   ACC_W   = 48
//
// Notes:
//   - During B/S load phase, load_en is asserted and S rows are loaded
//     into the seven cores.
//   - During A/T compute phase, load_en is low and T rows are streamed
//     through the seven cores.
//   - All seven cores have identical latency, so Q1..Q7 remain aligned.
//=============================================================================

(* keep_hierarchy = "yes" *)
module smm1_cores #(
    parameter int N     = 8,     // BLOCK_N
    parameter int MW    = 17,    // T/S width
    parameter int ACC_W = 48     // Q width
) (
    input  logic                        clk,
    input  logic                        rst_n,

    //-------------------------------------------------------------------------
    // B/S load interface
    //-------------------------------------------------------------------------
    input  logic                        load_en,
    input  logic [$clog2(N)-1:0]        load_row,

    input  logic signed [N*MW-1:0]      s1,
    input  logic signed [N*MW-1:0]      s2,
    input  logic signed [N*MW-1:0]      s3,
    input  logic signed [N*MW-1:0]      s4,
    input  logic signed [N*MW-1:0]      s5,
    input  logic signed [N*MW-1:0]      s6,
    input  logic signed [N*MW-1:0]      s7,

    //-------------------------------------------------------------------------
    // A/T streaming interface
    //-------------------------------------------------------------------------
    input  logic                        t_valid,

    input  logic signed [N*MW-1:0]      t1,
    input  logic signed [N*MW-1:0]      t2,
    input  logic signed [N*MW-1:0]      t3,
    input  logic signed [N*MW-1:0]      t4,
    input  logic signed [N*MW-1:0]      t5,
    input  logic signed [N*MW-1:0]      t6,
    input  logic signed [N*MW-1:0]      t7,

    //-------------------------------------------------------------------------
    // Q output interface
    //-------------------------------------------------------------------------
    output logic                        q_valid,

    output logic signed [N*ACC_W-1:0]   q1,
    output logic signed [N*ACC_W-1:0]   q2,
    output logic signed [N*ACC_W-1:0]   q3,
    output logic signed [N*ACC_W-1:0]   q4,
    output logic signed [N*ACC_W-1:0]   q5,
    output logic signed [N*ACC_W-1:0]   q6,
    output logic signed [N*ACC_W-1:0]   q7
);

    //-------------------------------------------------------------------------
    // Internal valid signals from each core
    //-------------------------------------------------------------------------
    logic q1_valid;
    logic q2_valid;
    logic q3_valid;
    logic q4_valid;
    logic q5_valid;
    logic q6_valid;
    logic q7_valid;

    //-------------------------------------------------------------------------
    // Core 1: Q1 = T1 * S1
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q1 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s1),
        .a_valid(t_valid), .a_row(t1),
        .c_valid(q1_valid), .c_row(q1)
    );

    //-------------------------------------------------------------------------
    // Core 2: Q2 = T2 * S2
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q2 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s2),
        .a_valid(t_valid), .a_row(t2),
        .c_valid(q2_valid), .c_row(q2)
    );

    //-------------------------------------------------------------------------
    // Core 3: Q3 = T3 * S3
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q3 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s3),
        .a_valid(t_valid), .a_row(t3),
        .c_valid(q3_valid), .c_row(q3)
    );

    //-------------------------------------------------------------------------
    // Core 4: Q4 = T4 * S4
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q4 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s4),
        .a_valid(t_valid), .a_row(t4),
        .c_valid(q4_valid), .c_row(q4)
    );

    //-------------------------------------------------------------------------
    // Core 5: Q5 = T5 * S5
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q5 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s5),
        .a_valid(t_valid), .a_row(t5),
        .c_valid(q5_valid), .c_row(q5)
    );

    //-------------------------------------------------------------------------
    // Core 6: Q6 = T6 * S6
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q6 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s6),
        .a_valid(t_valid), .a_row(t6),
        .c_valid(q6_valid), .c_row(q6)
    );

    //-------------------------------------------------------------------------
    // Core 7: Q7 = T7 * S7
    //-------------------------------------------------------------------------
    mm_ws_systolic #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q7 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s7),
        .a_valid(t_valid), .a_row(t7),
        .c_valid(q7_valid), .c_row(q7)
    );

    //-------------------------------------------------------------------------
    // All seven cores have the same pipeline latency,
    // therefore Q1..Q7 are aligned.
    //-------------------------------------------------------------------------
    assign q_valid = q1_valid & q2_valid & q3_valid & q4_valid
                   & q5_valid & q6_valid & q7_valid;

endmodule
