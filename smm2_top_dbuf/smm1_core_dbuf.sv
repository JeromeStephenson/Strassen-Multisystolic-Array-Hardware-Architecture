//=============================================================================
// smm1_core_dbuf.sv
// Reusable one-level Strassen datapath core with Double-Buffering support.
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

    input  logic                                   b_valid,
    input  logic [$clog2(BASE_N)-1:0]              b_row_idx,
    input  logic signed [4*BASE_N*IN_W-1:0]        b_row,

    input  logic                                   a_valid,
    input  logic signed [4*BASE_N*IN_W-1:0]        a_row,

    // Double-buffer control from top-level FSM
    input  logic                                   swap,
    output logic                                   active_bank,

    output logic                                   core_load_en,
    output logic [$clog2(BASE_N)-1:0]              core_load_row,

    output logic                                   c_valid,
    output logic signed [4*BASE_N*OUT_W-1:0]       c_row
);

    logic t_valid, s_valid, q_valid, c_valid_post;
    logic signed [BASE_N*MID_W-1:0] t1, t2, t3, t4, t5, t6, t7;
    logic signed [BASE_N*MID_W-1:0] s1, s2, s3, s4, s5, s6, s7;
    logic signed [BASE_N*ACC_W-1:0] q1, q2, q3, q4, q5, q6, q7;
    logic signed [BASE_N*OUT_W-1:0] c11_row, c12_row, c21_row, c22_row;
    logic [$clog2(BASE_N)-1:0] load_row_d;

    //-------------------------------------------------------------------------
    // Double Buffer Bank Tracking
    // Toggles exactly when the top-level FSM asserts the 1-cycle 'swap' pulse
    //-------------------------------------------------------------------------
    logic current_bank;
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) current_bank <= 1'b0;
        else if (swap) current_bank <= ~current_bank;
    end
    assign active_bank = current_bank;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) load_row_d <= '0;
        else if (b_valid) load_row_d <= b_row_idx;
    end

    assign core_load_en  = s_valid;
    assign core_load_row = load_row_d;

    strassen_pre #(.G(1), .N(BASE_N), .W(IN_W), .MW(MID_W)) u_strassen_pre (
        .clk(clk), .rst_n(rst_n), .a_valid(a_valid), .b_valid(b_valid),
        .a_row(a_row), .b_row(b_row),
        .t1(t1), .t2(t2), .t3(t3), .t4(t4), .t5(t5), .t6(t6), .t7(t7),
        .s1(s1), .s2(s2), .s3(s3), .s4(s4), .s5(s5), .s6(s6), .s7(s7),
        .t_valid(t_valid), .s_valid(s_valid)
    );

    smm1_cores_dbuf #(.N(BASE_N), .MW(MID_W), .ACC_W(ACC_W)) u_smm1_cores_dbuf (
        .clk(clk), .rst_n(rst_n),
        .a_bank(current_bank), .b_bank(~current_bank), // A reads active, B writes to shadow
        .load_en(core_load_en), .load_row(core_load_row),
        .s1(s1), .s2(s2), .s3(s3), .s4(s4), .s5(s5), .s6(s6), .s7(s7),
        .t_valid(t_valid),
        .t1(t1), .t2(t2), .t3(t3), .t4(t4), .t5(t5), .t6(t6), .t7(t7),
        .q_valid(q_valid),
        .q1(q1), .q2(q2), .q3(q3), .q4(q4), .q5(q5), .q6(q6), .q7(q7)
    );

    strassen_post #(.G(1), .N(BASE_N), .QW(ACC_W), .CW(OUT_W)) u_strassen_post (
        .clk(clk), .rst_n(rst_n), .in_valid(q_valid),
        .q1(q1), .q2(q2), .q3(q3), .q4(q4), .q5(q5), .q6(q6), .q7(q7),
        .c11(c11_row), .c12(c12_row), .c21(c21_row), .c22(c22_row), .out_valid(c_valid_post)
    );

    assign c_valid = c_valid_post;
    assign c_row   = {c22_row, c21_row, c12_row, c11_row};

endmodule
