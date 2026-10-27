//=============================================================================
// mm_ws_systolic_dbuf.sv
//
// Double-buffered weight-stationary pipelined systolic matrix multiplication
// core.
//
// Function:
//   C = A * B
//
// Paper-style feature (extra b buffer):
//   - Active b buffer is used for the current matrix multiplication.
//   - Shadow b buffer can be loaded with the next B tile while the active
//     buffer is being used.
//   - A one-cycle swap pulse exchanges active and shadow buffers.
//
// For Strassen SMM1:
//   A-side input is one of T1..T7.
//   B-side loaded matrix is one of S1..S7.
//   Output is one of Q1..Q7.
//
// Target parameters:
//   X = 8, Y = 8, A_W = 17, B_W = 17, ACC_W = 48
//
// Datapath latency from a_valid to c_valid: X + 1 cycles
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm_ws_systolic_dbuf #(
    parameter int X      = 8,     // inner dimension / number of B rows
    parameter int Y      = 8,     // output width / number of B columns
    parameter int A_W    = 17,    // A/T input bit width
    parameter int B_W    = 17,    // B/S input bit width
    parameter int ACC_W  = 48     // accumulator/output bit width
) (
    input  logic                        clk,
    input  logic                        rst_n,

    //-------------------------------------------------------------------------
    // B/S shadow-buffer load interface
    //-------------------------------------------------------------------------
    input  logic                        load_en,
    input  logic [$clog2(X)-1:0]        load_row,
    input  logic signed [Y*B_W-1:0]     b_row,

    //-------------------------------------------------------------------------
    // Buffer swap control
    //
    // swap        : one-cycle pulse, exchanges active and shadow buffers
    // active_bank : 0 = b0 active, 1 = b1 active
    //-------------------------------------------------------------------------
    input  logic                        swap,
    output logic                        active_bank,

    //-------------------------------------------------------------------------
    // A/T streaming interface
    //-------------------------------------------------------------------------
    input  logic                        a_valid,
    input  logic signed [X*A_W-1:0]     a_row,

    //-------------------------------------------------------------------------
    // C/Q output interface
    //-------------------------------------------------------------------------
    output logic                        c_valid,
    output logic signed [Y*ACC_W-1:0]   c_row
);

    //-------------------------------------------------------------------------
    // Active bank register
    //-------------------------------------------------------------------------
    logic active_q;

    assign active_bank = active_q;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            active_q <= 1'b0;
        end else begin
            if (swap) begin
                active_q <= ~active_q;
            end
        end
    end

    //-------------------------------------------------------------------------
    // Double B/S storage
    //
    //   active_q == 0 : b0 active, b1 shadow
    //   active_q == 1 : b1 active, b0 shadow
    //-------------------------------------------------------------------------
    logic signed [B_W-1:0] b0 [X][Y];
    logic signed [B_W-1:0] b1 [X][Y];

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int k = 0; k < X; k++)
                for (int j = 0; j < Y; j++) begin
                    b0[k][j] <= '0;
                    b1[k][j] <= '0;
                end
        end else begin
            if (load_en) begin
                // Load into shadow bank only.
                if (active_q == 1'b0) begin
                    for (int j = 0; j < Y; j++)
                        b1[load_row][j] <= $signed(b_row[j*B_W +: B_W]);
                end else begin
                    for (int j = 0; j < Y; j++)
                        b0[load_row][j] <= $signed(b_row[j*B_W +: B_W]);
                end
            end
        end
    end

    //-------------------------------------------------------------------------
    // Active B/S mux for the MAC pipeline
    //
    // Important: swap must only be asserted when the current pipeline is
    // safe to switch (after the current tile has drained).
    //-------------------------------------------------------------------------
    logic signed [B_W-1:0] b_sel [X][Y];

    always_comb begin
        for (int k = 0; k < X; k++)
            for (int j = 0; j < Y; j++)
                b_sel[k][j] = active_q ? b1[k][j] : b0[k][j];
    end

    //-------------------------------------------------------------------------
    // Unpack A/T input row
    //-------------------------------------------------------------------------
    logic signed [A_W-1:0] a_in [X];

    always_comb begin
        for (int k = 0; k < X; k++)
            a_in[k] = $signed(a_row[k*A_W +: A_W]);
    end

    //-------------------------------------------------------------------------
    // A/T skewing network
    //-------------------------------------------------------------------------
    logic signed [A_W-1:0] a_skew [X][X+1];

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int k = 0; k < X; k++)
                for (int d = 0; d <= X; d++)
                    a_skew[k][d] <= '0;
        end else begin
            for (int k = 0; k < X; k++) begin
                a_skew[k][0] <= a_valid ? a_in[k] : '0;
                for (int d = 1; d <= X; d++)
                    a_skew[k][d] <= a_skew[k][d-1];
            end
        end
    end

    //-------------------------------------------------------------------------
    // Valid pipeline
    //-------------------------------------------------------------------------
    logic [X:0] v_skew;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            v_skew <= '0;
        end else begin
            v_skew[0] <= a_valid;
            for (int d = 1; d <= X; d++)
                v_skew[d] <= v_skew[d-1];
        end
    end

    //-------------------------------------------------------------------------
    // Weight-stationary partial-sum pipeline
    //
    //   p[0][j]   = 0
    //   p[k+1][j] = p[k][j] + A[k] * B_sel[k][j]
    //   p[X][j]   = sum_k A[k]*B_sel[k][j]
    //-------------------------------------------------------------------------
    logic signed [ACC_W-1:0] p [X+1][Y];

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int s = 0; s <= X; s++)
                for (int j = 0; j < Y; j++)
                    p[s][j] <= '0;
        end else begin
            for (int j = 0; j < Y; j++)
                p[0][j] <= '0;

            for (int k = 0; k < X; k++)
                for (int j = 0; j < Y; j++)
                    p[k+1][j] <= p[k][j]
                               + ($signed(a_skew[k][k]) * $signed(b_sel[k][j]));
        end
    end

    //-------------------------------------------------------------------------
    // Output
    //-------------------------------------------------------------------------
    assign c_valid = v_skew[X];

    for (genvar gj = 0; gj < Y; gj++) begin : g_c_out
        assign c_row[gj*ACC_W +: ACC_W] = p[X][gj];
    end

endmodule
