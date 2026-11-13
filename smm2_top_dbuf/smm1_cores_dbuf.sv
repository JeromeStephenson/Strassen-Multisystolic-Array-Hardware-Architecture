//=============================================================================
// smm1_cores_dbuf.sv
// Seven-core Strassen MM array with Double-Buffering support.
//=============================================================================
(* keep_hierarchy = "yes" *)
module smm1_cores_dbuf #(
    parameter int N     = 8,
    parameter int MW    = 17,
    parameter int ACC_W = 48
) (
    input  logic                        clk,
    input  logic                        rst_n,
    input  logic                        a_bank,
    input  logic                        b_bank,
    input  logic                        load_en,
    input  logic [$clog2(N)-1:0]        load_row,
    input  logic signed [N*MW-1:0]      s1, s2, s3, s4, s5, s6, s7,
    input  logic                        t_valid,
    input  logic signed [N*MW-1:0]      t1, t2, t3, t4, t5, t6, t7,
    output logic                        q_valid,
    output logic signed [N*ACC_W-1:0]   q1, q2, q3, q4, q5, q6, q7
);

    logic q1_valid, q2_valid, q3_valid, q4_valid, q5_valid, q6_valid, q7_valid;

    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q1 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s1), .a_valid(t_valid), .a_row(t1), .c_valid(q1_valid), .c_row(q1));
    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q2 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s2), .a_valid(t_valid), .a_row(t2), .c_valid(q2_valid), .c_row(q2));
    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q3 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s3), .a_valid(t_valid), .a_row(t3), .c_valid(q3_valid), .c_row(q3));
    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q4 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s4), .a_valid(t_valid), .a_row(t4), .c_valid(q4_valid), .c_row(q4));
    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q5 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s5), .a_valid(t_valid), .a_row(t5), .c_valid(q5_valid), .c_row(q5));
    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q6 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s6), .a_valid(t_valid), .a_row(t6), .c_valid(q6_valid), .c_row(q6));
    mm_ws_systolic_dbuf #(.X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)) u_q7 (.clk(clk), .rst_n(rst_n), .a_bank(a_bank), .b_bank(b_bank), .load_en(load_en), .load_row(load_row), .b_row(s7), .a_valid(t_valid), .a_row(t7), .c_valid(q7_valid), .c_row(q7));

    assign q_valid = q1_valid & q2_valid & q3_valid & q4_valid & q5_valid & q6_valid & q7_valid;

endmodule
