//=============================================================================
// mxu_synth_wrap.sv
//
// Post-implementation (place & route) measurement wrapper.
//
// The wide A/B/C buses are driven/sunk by INTERNAL registers so the design
// can be implemented without thousands of package pins. Only clk/rst_n/
// start/b_valid/a_valid/done are real I/O.
//
// DUT_SEL selects which MXU top to measure:
//   0 = mm1_top          4 = mm2_top_dbuf
//   1 = smm1_top         5 = smm2_top_dbuf
//   2 = mm2_top          6 = mm2_top_gemm
//   3 = smm2_top         7 = smm2_top_gemm
//
// To read results: use report_utilization -hierarchical and take the u_dut
// subtree for LUT/FF/DSP. Read Fmax from report_timing_summary.
//=============================================================================

(* keep_hierarchy = "yes" *)
module mxu_synth_wrap #(
    parameter int DUT_SEL   = 3,   // see table above
    parameter int NUM_TILES = 2    // used only by the GEMM tops
) (
    input  logic clk,
    input  logic rst_n,
    input  logic start,
    input  logic b_valid,
    input  logic a_valid,
    output logic done
);

    //-------------------------------------------------------------------------
    // Bus widths per design family
    //-------------------------------------------------------------------------
    localparam int IN_ROW_BITS  = (DUT_SEL <= 1) ? 512  : 2048;
    localparam int OUT_ROW_BITS = (DUT_SEL <= 1) ? 1600 : 6656;
    localparam int TILE_CNT_W   = 16;

    //-------------------------------------------------------------------------
    // Free-running counter -> non-constant data inputs.
    // Keeping the inputs non-constant prevents Vivado from constant-propagating
    // the multipliers/adders away, so all DSPs are preserved.
    //-------------------------------------------------------------------------
    (* DONT_TOUCH = "yes" *) logic [15:0] cnt;
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) cnt <= 16'd1;
        else        cnt <= cnt + 16'd1;
    end

    (* DONT_TOUCH = "yes" *) logic [IN_ROW_BITS-1:0] a_row_r;
    (* DONT_TOUCH = "yes" *) logic [IN_ROW_BITS-1:0] b_row_r;
    always_ff @(posedge clk) begin
        a_row_r <= { (IN_ROW_BITS/8) { cnt[7:0]  } };
        b_row_r <= { (IN_ROW_BITS/8) { cnt[15:8] } };
    end

    //-------------------------------------------------------------------------
    // Sink the wide C output so the post-addition network is not optimized away
    //-------------------------------------------------------------------------
    logic [OUT_ROW_BITS-1:0] c_row_w;
    (* DONT_TOUCH = "yes" *) logic [OUT_ROW_BITS-1:0] c_row_r;
    always_ff @(posedge clk) c_row_r <= c_row_w;

    //-------------------------------------------------------------------------
    // GEMM-only control
    //-------------------------------------------------------------------------
    logic active_bank_w;
    logic [TILE_CNT_W-1:0] num_tiles_w;
    assign num_tiles_w = TILE_CNT_W'(NUM_TILES);

    //-------------------------------------------------------------------------
    // DUT instantiation selected by DUT_SEL
    //-------------------------------------------------------------------------
    generate
        if (DUT_SEL == 0) begin : g_dut
            mm1_top u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .busy(), .done(done),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 1) begin : g_dut
            smm1_top u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .busy(), .done(done),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 2) begin : g_dut
            mm2_top u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .busy(), .done(done),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 3) begin : g_dut
            smm2_top u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .busy(), .done(done),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 4) begin : g_dut
            mm2_top_dbuf u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .busy(), .done(done), .active_bank(active_bank_w),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 5) begin : g_dut
            smm2_top_dbuf u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .busy(), .done(done), .active_bank(active_bank_w),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 6) begin : g_dut
            mm2_top_gemm u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .num_tiles(num_tiles_w),
                .busy(), .done(done), .active_bank(active_bank_w),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
        else if (DUT_SEL == 7) begin : g_dut
            smm2_top_gemm u_dut (
                .clk(clk), .rst_n(rst_n), .start(start),
                .num_tiles(num_tiles_w),
                .busy(), .done(done), .active_bank(active_bank_w),
                .b_ready(), .b_valid(b_valid), .b_row(b_row_r),
                .a_ready(), .a_valid(a_valid), .a_row(a_row_r),
                .c_valid(), .c_row(c_row_w)
            );
        end
    endgenerate

endmodule
