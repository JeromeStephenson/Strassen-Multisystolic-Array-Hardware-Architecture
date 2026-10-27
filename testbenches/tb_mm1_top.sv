//=============================================================================
// tb_mm1_top.sv
//
// Self-checking testbench for the baseline conventional multisystolic top
// (mm1_top.sv). Same tests as tb_smm1_top_dbuf so SMM1 and MM1 can be
// compared functionally and in synthesis.
//=============================================================================

`timescale 1ns / 1ps

module tb_mm1_top;

    //-------------------------------------------------------------------------
    // Parameters (must match DUT)
    //-------------------------------------------------------------------------
    localparam int TOP_N    = 16;
    localparam int SMM_R    = 1;
    localparam int BLOCK_N  = 8;
    localparam int DATA_W   = 16;
    localparam int ACC_W    = 48;
    localparam int OUT_W    = ACC_W + 2;            // 50

    localparam int IN_ROW_BITS   = 4 * BLOCK_N * DATA_W;   // 512
    localparam int OUT_ROW_BITS  = 4 * BLOCK_N * OUT_W;    // 1600
    localparam int QUAD_IN_BITS  = BLOCK_N * DATA_W;       // 128
    localparam int QUAD_OUT_BITS = BLOCK_N * OUT_W;        // 400

    //-------------------------------------------------------------------------
    // Clock
    //-------------------------------------------------------------------------
    logic clk = 1'b0;
    always #5 clk = ~clk;

    //-------------------------------------------------------------------------
    // DUT interface signals
    //-------------------------------------------------------------------------
    logic                       rst_n;
    logic                       start;
    logic                       busy;
    logic                       done;
    logic                       b_ready;
    logic                       b_valid;
    logic [IN_ROW_BITS-1:0]     b_row;
    logic                       a_ready;
    logic                       a_valid;
    logic [IN_ROW_BITS-1:0]     a_row;
    logic                       c_valid;
    logic [OUT_ROW_BITS-1:0]    c_row;

    //-------------------------------------------------------------------------
    // DUT (baseline MM1 - note: no active_bank port)
    //-------------------------------------------------------------------------
    mm1_top #(
        .TOP_N   (TOP_N),
        .SMM_R   (SMM_R),
        .BLOCK_N (BLOCK_N),
        .DATA_W  (DATA_W),
        .ACC_W   (ACC_W),
        .OUT_W   (OUT_W)
    ) dut (
        .clk         (clk),
        .rst_n       (rst_n),
        .start       (start),
        .busy        (busy),
        .done        (done),
        .b_ready     (b_ready),
        .b_valid     (b_valid),
        .b_row       (b_row),
        .a_ready     (a_ready),
        .a_valid     (a_valid),
        .a_row       (a_row),
        .c_valid     (c_valid),
        .c_row       (c_row)
    );

    //-------------------------------------------------------------------------
    // Matrices
    //-------------------------------------------------------------------------
    logic signed [DATA_W-1:0] A [TOP_N][TOP_N];
    logic signed [DATA_W-1:0] B [TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_expected [TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_dut      [TOP_N][TOP_N];

    int c_idx;
    int errors;

    //-------------------------------------------------------------------------
    // Capture C block-rows from DUT
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c_idx <= 0;
        end else if (c_valid && (c_idx < BLOCK_N)) begin
            for (int j = 0; j < BLOCK_N; j++) begin
                C_dut[c_idx][j]                 <= c_row[0*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
                C_dut[c_idx][BLOCK_N+j]         <= c_row[1*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
                C_dut[BLOCK_N+c_idx][j]         <= c_row[2*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
                C_dut[BLOCK_N+c_idx][BLOCK_N+j] <= c_row[3*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
            end
            c_idx <= c_idx + 1;
        end
    end

    //-------------------------------------------------------------------------
    // Watchdog
    //-------------------------------------------------------------------------
    initial begin
        #2_000_000;
        $display("TB ERROR: global watchdog timeout - simulation hung.");
        $finish;
    end

    //-------------------------------------------------------------------------
    // Task: reset
    //-------------------------------------------------------------------------
    task automatic reset_dut();
        begin
            rst_n   = 1'b0;
            start   = 1'b0;
            b_valid = 1'b0;
            a_valid = 1'b0;
            b_row   = '0;
            a_row   = '0;
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++)
                    C_dut[r][c] = '0;
            repeat (10) @(posedge clk);
            rst_n = 1'b1;
            repeat (5) @(posedge clk);
        end
    endtask

    //-------------------------------------------------------------------------
    // Task: golden model C = A * B
    //-------------------------------------------------------------------------
    task automatic compute_expected();
        longint sum;
        begin
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    sum = 0;
                    for (int k = 0; k < TOP_N; k++)
                        sum = sum + (longint'(A[r][k]) * longint'(B[k][c]));
                    C_expected[r][c] = sum;
                end
        end
    endtask

    //-------------------------------------------------------------------------
    // Test initializers
    //-------------------------------------------------------------------------
    task automatic init_identity_test();
        int val;
        begin
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++)
                    A[r][c] = (r == c) ? 16'sd1 : 16'sd0;
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    val     = (r * TOP_N + c) % 201;
                    B[r][c] = val - 100;
                end
            compute_expected();
        end
    endtask

    task automatic init_small_ramp_test();
        int val;
        begin
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    val     = ((r * TOP_N + c) % 17) - 8;
                    A[r][c] = val;
                    val     = ((r * TOP_N + c) % 13) - 6;
                    B[r][c] = val;
                end
            compute_expected();
        end
    endtask

    //-------------------------------------------------------------------------
    // Packing functions
    //-------------------------------------------------------------------------
    function automatic logic [IN_ROW_BITS-1:0] get_a_row(input int i);
        logic [QUAD_IN_BITS-1:0] a11, a12, a21, a22;
        begin
            a11 = '0; a12 = '0; a21 = '0; a22 = '0;
            for (int j = 0; j < BLOCK_N; j++) begin
                a11[j*DATA_W +: DATA_W] = A[i][j];
                a12[j*DATA_W +: DATA_W] = A[i][BLOCK_N+j];
                a21[j*DATA_W +: DATA_W] = A[BLOCK_N+i][j];
                a22[j*DATA_W +: DATA_W] = A[BLOCK_N+i][BLOCK_N+j];
            end
            get_a_row = {a22, a21, a12, a11};
        end
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_b_row(input int i);
        logic [QUAD_IN_BITS-1:0] b11, b12, b21, b22;
        begin
            b11 = '0; b12 = '0; b21 = '0; b22 = '0;
            for (int j = 0; j < BLOCK_N; j++) begin
                b11[j*DATA_W +: DATA_W] = B[i][j];
                b12[j*DATA_W +: DATA_W] = B[i][BLOCK_N+j];
                b21[j*DATA_W +: DATA_W] = B[BLOCK_N+i][j];
                b22[j*DATA_W +: DATA_W] = B[BLOCK_N+i][BLOCK_N+j];
            end
            get_b_row = {b22, b21, b12, b11};
        end
    endfunction

    //-------------------------------------------------------------------------
    // Handshake tasks
    //-------------------------------------------------------------------------
    task automatic send_b_row(input logic [IN_ROW_BITS-1:0] data);
        begin
            @(negedge clk);
            b_valid = 1'b1;
            b_row   = data;
            wait (b_ready);
            @(posedge clk);
            @(negedge clk);
            b_valid = 1'b0;
        end
    endtask

    task automatic send_a_row(input logic [IN_ROW_BITS-1:0] data);
        begin
            @(negedge clk);
            a_valid = 1'b1;
            a_row   = data;
            wait (a_ready);
            @(posedge clk);
            @(negedge clk);
            a_valid = 1'b0;
        end
    endtask

    //-------------------------------------------------------------------------
    // Wait for done
    //-------------------------------------------------------------------------
    task automatic wait_done();
        int t;
        begin
            t = 0;
            while (!done && t < 20000) begin
                @(posedge clk);
                t = t + 1;
            end
            if (!done) $fatal(1, "TB ERROR: timeout waiting for done.");
            repeat (10) @(posedge clk);
        end
    endtask

    //-------------------------------------------------------------------------
    // Compare
    //-------------------------------------------------------------------------
    task automatic compare_results(input string test_name);
        begin
            errors = 0;
            if (c_idx != BLOCK_N) begin
                $display("TB ERROR: %s captured %0d C rows, expected %0d.",
                         test_name, c_idx, BLOCK_N);
                errors = errors + 1;
            end
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++)
                    if (C_dut[r][c] !== C_expected[r][c]) begin
                        errors = errors + 1;
                        if (errors <= 20)
                            $display("TB ERROR: %s C[%0d][%0d] DUT=%0d EXP=%0d",
                                     test_name, r, c, C_dut[r][c], C_expected[r][c]);
                    end
            $display("==================================================");
            if (errors == 0) $display("%s : PASS", test_name);
            else             $display("%s : FAIL (%0d errors)", test_name, errors);
            $display("==================================================");
        end
    endtask

    //-------------------------------------------------------------------------
    // Run one test
    //-------------------------------------------------------------------------
    task automatic run_test(input string test_name);
        begin
            $display("TB: [%s] reset...", test_name);
            reset_dut();
            @(negedge clk); start = 1'b1;
            @(negedge clk); start = 1'b0;
            for (int i = 0; i < BLOCK_N; i++) send_b_row(get_b_row(i));
            for (int i = 0; i < BLOCK_N; i++) send_a_row(get_a_row(i));
            wait_done();
            compare_results(test_name);
        end
    endtask

    //-------------------------------------------------------------------------
    // Main
    //-------------------------------------------------------------------------
    initial begin
        rst_n   = 1'b0;
        start   = 1'b0;
        b_valid = 1'b0;
        a_valid = 1'b0;
        b_row   = '0;
        a_row   = '0;

        init_identity_test();
        run_test("MM1 Test 1: Identity A");

        init_small_ramp_test();
        run_test("MM1 Test 2: Small ramp");

        $display("TB: simulation finished.");
        $finish;
    end

endmodule
