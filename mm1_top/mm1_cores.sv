//=============================================================================
// mm1_cores.sv
//
// Eight-core conventional blocked matrix multiplication array (baseline MM1,
// paper Fig. 5 / Eq. 2).
//
//   Q1 = A11*B11   Q2 = A11*B12   Q3 = A12*B21   Q4 = A12*B22
//   Q5 = A21*B11   Q6 = A21*B12   Q7 = A22*B21   Q8 = A22*B22
//
// Note: baseline uses 16-bit multiplier inputs (no Strassen bit growth).
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm1_cores #(
    parameter int N     = 8,     // BLOCK_N
    parameter int W     = 16,    // DATA_W (no growth for baseline)
    parameter int ACC_W = 48
) (
    input  logic                        clk,
    input  logic                        rst_n,

    input  logic                        load_en,
    input  logic [$clog2(N)-1:0]        load_row,

    // A quadrant rows
    input  logic signed [N*W-1:0]       a11,
    input  logic signed [N*W-1:0]       a12,
    input  logic signed [N*W-1:0]       a21,
    input  logic signed [N*W-1:0]       a22,

    // B quadrant rows
    input  logic signed [N*W-1:0]       b11,
    input  logic signed [N*W-1:0]       b12,
    input  logic signed [N*W-1:0]       b21,
    input  logic signed [N*W-1:0]       b22,

    input  logic                        t_valid,

    output logic                        q_valid,
    output logic signed [N*ACC_W-1:0]   q1,
    output logic signed [N*ACC_W-1:0]   q2,
    output logic signed [N*ACC_W-1:0]   q3,
    output logic signed [N*ACC_W-1:0]   q4,
    output logic signed [N*ACC_W-1:0]   q5,
    output logic signed [N*ACC_W-1:0]   q6,
    output logic signed [N*ACC_W-1:0]   q7,
    output logic signed [N*ACC_W-1:0]   q8
);

    logic q1_valid, q2_valid, q3_valid, q4_valid;
    logic q5_valid, q6_valid, q7_valid, q8_valid;

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q1 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b11), .a_valid(t_valid), .a_row(a11),
        .c_valid(q1_valid), .c_row(q1));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q2 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b12), .a_valid(t_valid), .a_row(a11),
        .c_valid(q2_valid), .c_row(q2));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q3 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b21), .a_valid(t_valid), .a_row(a12),
        .c_valid(q3_valid), .c_row(q3));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q4 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b22), .a_valid(t_valid), .a_row(a12),
        .c_valid(q4_valid), .c_row(q4));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q5 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b11), .a_valid(t_valid), .a_row(a21),
        .c_valid(q5_valid), .c_row(q5));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q6 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b12), .a_valid(t_valid), .a_row(a21),
        .c_valid(q6_valid), .c_row(q6));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q7 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b21), .a_valid(t_valid), .a_row(a22),
        .c_valid(q7_valid), .c_row(q7));

    mm_ws_systolic #(.X(N), .Y(N), .A_W(W), .B_W(W), .ACC_W(ACC_W)) u_q8 (
        .clk(clk), .rst_n(rst_n), .load_en(load_en), .load_row(load_row),
        .b_row(b22), .a_valid(t_valid), .a_row(a22),
        .c_valid(q8_valid), .c_row(q8));

    assign q_valid = q1_valid & q2_valid & q3_valid & q4_valid
                   & q5_valid & q6_valid & q7_valid & q8_valid;

endmodule
