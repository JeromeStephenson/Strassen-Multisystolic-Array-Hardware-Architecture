
//=============================================================================
// strassen_pre_generic.sv
//
// Generic Strassen pre-addition network.
//
// This module implements the Strassen T/S formation for one recursion level:
//
//   T1 = A11 + A22       S1 = B11 + B22
//   T2 = A21 + A22       S2 = B11
//   T3 = A11             S3 = B12 - B22
//   T4 = A22             S4 = B21 - B11
//   T5 = A11 + A12       S5 = B22
//   T6 = A21 - A11       S6 = B11 + B12
//   T7 = A12 - A22       S7 = B21 + B22
//
// The input is a packed row vector containing four quadrant vectors:
//
//   a_row = {A22 vector, A21 vector, A12 vector, A11 vector}
//   b_row = {B22 vector, B21 vector, B12 vector, B11 vector}
//
// Each quadrant vector has G*N elements.
//
// For SMM1 top:
//   G = 1, N = 8
//
// For SMM2 top:
//   G = 4, N = 8
//
// Pipeline latency: 1 cycle.
//=============================================================================

(* keep_hierarchy = "yes" *)
module strassen_pre #(
    parameter int G  = 1,          // number of base-row groups per quadrant
    parameter int N  = 8,          // base row length
    parameter int W  = 16,         // input bit width
    parameter int MW = W + 1       // output bit width after Strassen pre-add
) (
    input  logic                              clk,
    input  logic                              rst_n,

    input  logic                              a_valid,
    input  logic                              b_valid,

    input  logic signed [4*G*N*W-1:0]         a_row,
    input  logic signed [4*G*N*W-1:0]         b_row,

    output logic signed [G*N*MW-1:0]          t1,
    output logic signed [G*N*MW-1:0]          t2,
    output logic signed [G*N*MW-1:0]          t3,
    output logic signed [G*N*MW-1:0]          t4,
    output logic signed [G*N*MW-1:0]          t5,
    output logic signed [G*N*MW-1:0]          t6,
    output logic signed [G*N*MW-1:0]          t7,

    output logic signed [G*N*MW-1:0]          s1,
    output logic signed [G*N*MW-1:0]          s2,
    output logic signed [G*N*MW-1:0]          s3,
    output logic signed [G*N*MW-1:0]          s4,
    output logic signed [G*N*MW-1:0]          s5,
    output logic signed [G*N*MW-1:0]          s6,
    output logic signed [G*N*MW-1:0]          s7,

    output logic                              t_valid,
    output logic                              s_valid
);

    //-------------------------------------------------------------------------
    // Local parameters
    //-------------------------------------------------------------------------
    localparam int M  = G * N;         // elements per quadrant vector
    localparam int QW = M * W;         // quadrant input width in bits
    localparam int OW = M * MW;        // quadrant output width in bits

    //-------------------------------------------------------------------------
    // Unpacked quadrant elements
    //-------------------------------------------------------------------------
    logic signed [W-1:0] A11 [M];
    logic signed [W-1:0] A12 [M];
    logic signed [W-1:0] A21 [M];
    logic signed [W-1:0] A22 [M];

    logic signed [W-1:0] B11 [M];
    logic signed [W-1:0] B12 [M];
    logic signed [W-1:0] B21 [M];
    logic signed [W-1:0] B22 [M];

    //-------------------------------------------------------------------------
    // Valid-gated quadrant elements
    //-------------------------------------------------------------------------
    logic signed [W-1:0] A11g [M];
    logic signed [W-1:0] A12g [M];
    logic signed [W-1:0] A21g [M];
    logic signed [W-1:0] A22g [M];

    logic signed [W-1:0] B11g [M];
    logic signed [W-1:0] B12g [M];
    logic signed [W-1:0] B21g [M];
    logic signed [W-1:0] B22g [M];

    //-------------------------------------------------------------------------
    // Unpack A quadrants
    //
    // Packing order:
    //   a_row = {A22, A21, A12, A11}
    //-------------------------------------------------------------------------
    always_comb begin
        for (int i = 0; i < M; i++) begin
            A11[i] = $signed(a_row[0*QW + i*W +: W]);
            A12[i] = $signed(a_row[1*QW + i*W +: W]);
            A21[i] = $signed(a_row[2*QW + i*W +: W]);
            A22[i] = $signed(a_row[3*QW + i*W +: W]);
        end
    end

    //-------------------------------------------------------------------------
    // Unpack B quadrants
    //
    // Packing order:
    //   b_row = {B22, B21, B12, B11}
    //-------------------------------------------------------------------------
    always_comb begin
        for (int i = 0; i < M; i++) begin
            B11[i] = $signed(b_row[0*QW + i*W +: W]);
            B12[i] = $signed(b_row[1*QW + i*W +: W]);
            B21[i] = $signed(b_row[2*QW + i*W +: W]);
            B22[i] = $signed(b_row[3*QW + i*W +: W]);
        end
    end

    //-------------------------------------------------------------------------
    // Gate inputs with valid signals
    //-------------------------------------------------------------------------
    always_comb begin
        for (int i = 0; i < M; i++) begin
            A11g[i] = a_valid ? A11[i] : '0;
            A12g[i] = a_valid ? A12[i] : '0;
            A21g[i] = a_valid ? A21[i] : '0;
            A22g[i] = a_valid ? A22[i] : '0;

            B11g[i] = b_valid ? B11[i] : '0;
            B12g[i] = b_valid ? B12[i] : '0;
            B21g[i] = b_valid ? B21[i] : '0;
            B22g[i] = b_valid ? B22[i] : '0;
        end
    end

    //-------------------------------------------------------------------------
    // Sign extension W -> MW bits
    //-------------------------------------------------------------------------
    function automatic logic signed [MW-1:0] sext(
        input logic signed [W-1:0] x
    );
        if (MW == W) begin
            sext = x;
        end else begin
            sext = $signed({{(MW-W){x[W-1]}}, x});
        end
    endfunction

    //-------------------------------------------------------------------------
    // Strassen pre-addition pipeline
    //
    // One registered stage.
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            t_valid <= 1'b0;
            s_valid <= 1'b0;

            t1 <= '0;
            t2 <= '0;
            t3 <= '0;
            t4 <= '0;
            t5 <= '0;
            t6 <= '0;
            t7 <= '0;

            s1 <= '0;
            s2 <= '0;
            s3 <= '0;
            s4 <= '0;
            s5 <= '0;
            s6 <= '0;
            s7 <= '0;
        end else begin
            t_valid <= a_valid;
            s_valid <= b_valid;

            for (int i = 0; i < M; i++) begin
                //-------------------------------------------------------------
                // T matrices
                //-------------------------------------------------------------
                t1[i*MW +: MW] <= sext(A11g[i]) + sext(A22g[i]);
                t2[i*MW +: MW] <= sext(A21g[i]) + sext(A22g[i]);
                t3[i*MW +: MW] <= sext(A11g[i]);
                t4[i*MW +: MW] <= sext(A22g[i]);
                t5[i*MW +: MW] <= sext(A11g[i]) + sext(A12g[i]);
                t6[i*MW +: MW] <= sext(A21g[i]) - sext(A11g[i]);
                t7[i*MW +: MW] <= sext(A12g[i]) - sext(A22g[i]);

                //-------------------------------------------------------------
                // S matrices
                //-------------------------------------------------------------
                s1[i*MW +: MW] <= sext(B11g[i]) + sext(B22g[i]);
                s2[i*MW +: MW] <= sext(B11g[i]);
                s3[i*MW +: MW] <= sext(B12g[i]) - sext(B22g[i]);
                s4[i*MW +: MW] <= sext(B21g[i]) - sext(B11g[i]);
                s5[i*MW +: MW] <= sext(B22g[i]);
                s6[i*MW +: MW] <= sext(B11g[i]) + sext(B12g[i]);
                s7[i*MW +: MW] <= sext(B21g[i]) + sext(B22g[i]);
            end
        end
    end

endmodule
