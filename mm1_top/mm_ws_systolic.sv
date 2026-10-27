//=============================================================================
// mm_ws_systolic.sv
//
// Weight-stationary pipelined systolic matrix multiplication core.
//
// Function:
//   C = A * B
//
//   A is streamed one row per cycle (a_valid / a_row).
//   B is loaded row-by-row before computation (load_en / load_row / b_row).
//   c_valid is asserted X+1 cycles after a_valid.
//
// Default size for the 16x16 project: X = Y = 8.
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm_ws_systolic #(
    parameter int X      = 8,     // inner dimension / number of B rows
    parameter int Y      = 8,     // output width / number of B columns
    parameter int A_W    = 17,    // A/T input bit width
    parameter int B_W    = 17,    // B/S input bit width
    parameter int ACC_W  = 48     // accumulator/output bit width
) (
    input  logic                        clk,
    input  logic                        rst_n,

    //-------------------------------------------------------------------------
    // B/S load interface
    //-------------------------------------------------------------------------
    input  logic                        load_en,
    input  logic [$clog2(X)-1:0]        load_row,
    input  logic signed [Y*B_W-1:0]     b_row,

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
    // B/S storage inside PEs: b[k][j] stores B[k][j]
    //-------------------------------------------------------------------------
    logic signed [B_W-1:0] b [X][Y];

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int k = 0; k < X; k++)
                for (int j = 0; j < Y; j++)
                    b[k][j] <= '0;
        end else begin
            if (load_en) begin
                for (int j = 0; j < Y; j++)
                    b[load_row][j] <= $signed(b_row[j*B_W +: B_W]);
            end
        end
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
    //
    // a_skew[k][k] provides A[k] delayed by k+1 cycles so that it meets the
    // partial sum of row k at the right time inside the column pipeline.
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
    // Valid pipeline (v_skew[X] aligns with the final output row)
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
    //   p[k+1][j] = p[k][j] + A[k] * B[k][j]
    //   p[X][j]   = sum_k A[k]*B[k][j]
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
                               + ($signed(a_skew[k][k]) * $signed(b[k][j]));
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


//=============================================================================
// mm_ws_systolic.sv
//
// Weight-stationary pipelined systolic matrix multiplication core.
//
// Function:
//   C = A * B
//
//   A is streamed one row per cycle (a_valid / a_row).
//   B is loaded row-by-row before computation (load_en / load_row / b_row).
//   c_valid is asserted X+1 cycles after a_valid.
//
// Default size for the 16x16 project: X = Y = 8.
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm_ws_systolic #(
    parameter int X      = 8,     // inner dimension / number of B rows
    parameter int Y      = 8,     // output width / number of B columns
    parameter int A_W    = 17,    // A/T input bit width
    parameter int B_W    = 17,    // B/S input bit width
    parameter int ACC_W  = 48     // accumulator/output bit width
) (
    input  logic                        clk,
    input  logic                        rst_n,

    //-------------------------------------------------------------------------
    // B/S load interface
    //-------------------------------------------------------------------------
    input  logic                        load_en,
    input  logic [$clog2(X)-1:0]        load_row,
    input  logic signed [Y*B_W-1:0]     b_row,

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
    // B/S storage inside PEs: b[k][j] stores B[k][j]
    //-------------------------------------------------------------------------
    logic signed [B_W-1:0] b [X][Y];

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (int k = 0; k < X; k++)
                for (int j = 0; j < Y; j++)
                    b[k][j] <= '0;
        end else begin
            if (load_en) begin
                for (int j = 0; j < Y; j++)
                    b[load_row][j] <= $signed(b_row[j*B_W +: B_W]);
            end
        end
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
    //
    // a_skew[k][k] provides A[k] delayed by k+1 cycles so that it meets the
    // partial sum of row k at the right time inside the column pipeline.
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
    // Valid pipeline (v_skew[X] aligns with the final output row)
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
    //   p[k+1][j] = p[k][j] + A[k] * B[k][j]
    //   p[X][j]   = sum_k A[k]*B[k][j]
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
                               + ($signed(a_skew[k][k]) * $signed(b[k][j]));
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
