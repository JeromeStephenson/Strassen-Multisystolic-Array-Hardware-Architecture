//=============================================================================
// strassen_post_generic.sv
//
// Generic Strassen post-addition network.
//
// Implements:
//
//   C11 = Q1 + Q4 - Q5 + Q7
//   C12 = Q3 + Q5
//   C21 = Q2 + Q4
//   C22 = Q1 - Q2 + Q3 + Q6
//
// Each Q input is a packed vector of G*N elements.
//
// Pipeline: two stages.
//
// Stage 1:
//   C11_a = Q1 + Q4       C11_b = Q7 - Q5
//   C12_a = Q3 + Q5       C12_b = 0
//   C21_a = Q2 + Q4       C21_b = 0
//   C22_a = Q1 - Q2       C22_b = Q3 + Q6
//
// Stage 2:
//   C11 = C11_a + C11_b
//   C12 = C12_a + C12_b
//   C21 = C21_a + C21_b
//   C22 = C22_a + C22_b
//=============================================================================

(* keep_hierarchy = "yes" *)
module strassen_post #(
    parameter int G  = 1,
    parameter int N  = 8,
    parameter int QW = 48,
    parameter int CW = QW + 2
) (
    input  logic                              clk,
    input  logic                              rst_n,

    input  logic                              in_valid,

    input  logic signed [G*N*QW-1:0]          q1,
    input  logic signed [G*N*QW-1:0]          q2,
    input  logic signed [G*N*QW-1:0]          q3,
    input  logic signed [G*N*QW-1:0]          q4,
    input  logic signed [G*N*QW-1:0]          q5,
    input  logic signed [G*N*QW-1:0]          q6,
    input  logic signed [G*N*QW-1:0]          q7,

    output logic signed [G*N*CW-1:0]          c11,
    output logic signed [G*N*CW-1:0]          c12,
    output logic signed [G*N*CW-1:0]          c21,
    output logic signed [G*N*CW-1:0]          c22,

    output logic                              out_valid
);

    //-------------------------------------------------------------------------
    // Local parameters
    //-------------------------------------------------------------------------
    localparam int M  = G * N;
    localparam int IQW = M * QW;
    localparam int OCW = M * CW;

    //-------------------------------------------------------------------------
    // Unpacked Q elements
    //-------------------------------------------------------------------------
    logic signed [QW-1:0] Q1 [M];
    logic signed [QW-1:0] Q2 [M];
    logic signed [QW-1:0] Q3 [M];
    logic signed [QW-1:0] Q4 [M];
    logic signed [QW-1:0] Q5 [M];
    logic signed [QW-1:0] Q6 [M];
    logic signed [QW-1:0] Q7 [M];

    //-------------------------------------------------------------------------
    // Valid-gated Q elements
    //-------------------------------------------------------------------------
    logic signed [QW-1:0] Q1g [M];
    logic signed [QW-1:0] Q2g [M];
    logic signed [QW-1:0] Q3g [M];
    logic signed [QW-1:0] Q4g [M];
    logic signed [QW-1:0] Q5g [M];
    logic signed [QW-1:0] Q6g [M];
    logic signed [QW-1:0] Q7g [M];

    //-------------------------------------------------------------------------
    // Unpack Q vectors
    //-------------------------------------------------------------------------
    always_comb begin
        for (int i = 0; i < M; i++) begin
            Q1[i] = $signed(q1[i*QW +: QW]);
            Q2[i] = $signed(q2[i*QW +: QW]);
            Q3[i] = $signed(q3[i*QW +: QW]);
            Q4[i] = $signed(q4[i*QW +: QW]);
            Q5[i] = $signed(q5[i*QW +: QW]);
            Q6[i] = $signed(q6[i*QW +: QW]);
            Q7[i] = $signed(q7[i*QW +: QW]);
        end
    end

    //-------------------------------------------------------------------------
    // Gate with valid
    //-------------------------------------------------------------------------
    always_comb begin
        for (int i = 0; i < M; i++) begin
            Q1g[i] = in_valid ? Q1[i] : '0;
            Q2g[i] = in_valid ? Q2[i] : '0;
            Q3g[i] = in_valid ? Q3[i] : '0;
            Q4g[i] = in_valid ? Q4[i] : '0;
            Q5g[i] = in_valid ? Q5[i] : '0;
            Q6g[i] = in_valid ? Q6[i] : '0;
            Q7g[i] = in_valid ? Q7[i] : '0;
        end
    end

    //-------------------------------------------------------------------------
    // Sign extension QW -> CW
    //-------------------------------------------------------------------------
    function automatic logic signed [CW-1:0] sextq(
        input logic signed [QW-1:0] x
    );
        if (CW == QW) begin
            sextq = x;
        end else begin
            sextq = $signed({{(CW-QW){x[QW-1]}}, x});
        end
    endfunction

    //-------------------------------------------------------------------------
    // Stage 1 registers
    //-------------------------------------------------------------------------
    logic signed [OCW-1:0] c11_a;
    logic signed [OCW-1:0] c11_b;

    logic signed [OCW-1:0] c12_a;
    logic signed [OCW-1:0] c12_b;

    logic signed [OCW-1:0] c21_a;
    logic signed [OCW-1:0] c21_b;

    logic signed [OCW-1:0] c22_a;
    logic signed [OCW-1:0] c22_b;

    logic stage1_valid;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            stage1_valid <= 1'b0;

            c11_a <= '0;
            c11_b <= '0;

            c12_a <= '0;
            c12_b <= '0;

            c21_a <= '0;
            c21_b <= '0;

            c22_a <= '0;
            c22_b <= '0;
        end else begin
            stage1_valid <= in_valid;

            for (int i = 0; i < M; i++) begin
                c11_a[i*CW +: CW] <= sextq(Q1g[i]) + sextq(Q4g[i]);
                c11_b[i*CW +: CW] <= sextq(Q7g[i]) - sextq(Q5g[i]);

                c12_a[i*CW +: CW] <= sextq(Q3g[i]) + sextq(Q5g[i]);
                c12_b[i*CW +: CW] <= '0;

                c21_a[i*CW +: CW] <= sextq(Q2g[i]) + sextq(Q4g[i]);
                c21_b[i*CW +: CW] <= '0;

                c22_a[i*CW +: CW] <= sextq(Q1g[i]) - sextq(Q2g[i]);
                c22_b[i*CW +: CW] <= sextq(Q3g[i]) + sextq(Q6g[i]);
            end
        end
    end

    //-------------------------------------------------------------------------
    // Stage 2 registers
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            out_valid <= 1'b0;

            c11 <= '0;
            c12 <= '0;
            c21 <= '0;
            c22 <= '0;
        end else begin
            out_valid <= stage1_valid;

            for (int i = 0; i < M; i++) begin
                c11[i*CW +: CW] <= c11_a[i*CW +: CW] + c11_b[i*CW +: CW];
                c12[i*CW +: CW] <= c12_a[i*CW +: CW] + c12_b[i*CW +: CW];
                c21[i*CW +: CW] <= c21_a[i*CW +: CW] + c21_b[i*CW +: CW];
                c22[i*CW +: CW] <= c22_a[i*CW +: CW] + c22_b[i*CW +: CW];
            end
        end
    end


endmodule
