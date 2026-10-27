//=============================================================================
// tb_smm2_top_dbuf.sv
//
// Self-checking testbench for smm2_top_dbuf.sv (double-buffered SMM2).
//
// For a single tile the external handshake is identical to smm2_top:
//   1. Feed 8 top-level B rows while b_ready is high.
//   2. Feed 8 top-level A rows while a_ready is high.
//   3. Receive 8 top-level C rows on c_valid / c_row.
// The internal shadow-load + swap is handled by the DUT FSM.
//
// Tests:
//   1. A = identity, B = deterministic matrix -> C must equal B
//   2. A,B = small deterministic ramp matrices -> C = A*B
//=============================================================================

`timescale 1ns / 1ps

module tb_smm2_top_dbuf;

    //-------------------------------------------------------------------------
    // Parameters (must match DUT)
    //-------------------------------------------------------------------------
    localparam int TOP_N         = 32;
    localparam int SMM_R         = 2;
    localparam int BLOCK_N       = 8;
    localparam int TOP_G         = 4;

    localparam int DATA_W        = 16;
    localparam int TOP_MID_W     = DATA_W + 1;        // 17
    localparam int SUB_IN_W      = TOP_MID_W;         // 17
    localparam int SUB_MID_W     = SUB_IN_W + 1;      // 18
    localparam int ACC_W         = 48;
    localparam int SUB_OUT_W     = ACC_W + 2;         // 50
    localparam int OUT_W         = SUB_OUT_W + 2;     // 52

    //-------------------------------------------------------------------------
    // Packing parameters
    //-------------------------------------------------------------------------
    localparam int SUB_ELEMS       = BLOCK_N;                  // 8
    localparam int QUAD_ELEMS      = 4 * BLOCK_N;              // 32

    localparam int SUB_IN_BITS     = SUB_ELEMS * DATA_W;       // 128
    localparam int QUAD_IN_BITS    = QUAD_ELEMS * DATA_W;      // 512
    localparam int IN_ROW_BITS     = 4 * QUAD_IN_BITS;         // 2048

    localparam int SUB_OUT_BITS    = SUB_ELEMS * OUT_W;        // 416
    localparam int QUAD_OUT_BITS   = QUAD_ELEMS * OUT_W;       // 1664
    localparam int OUT_ROW_BITS    = 4 * QUAD_OUT_BITS;        // 6656

    //-------------------------------------------------------------------------
    // Clock
    //-------------------------------------------------------------------------
    logic clk = 1'b0;
    always #5 clk = ~clk;

    //-------------------------------------------------------------------------
    // DUT interface
    //-------------------------------------------------------------------------
    logic                            rst_n;
    logic                            start;
    logic                            busy;
    logic                            done;
    logic                            active_bank;

    logic                            b_ready;
    logic                            b_valid;
    logic signed [IN_ROW_BITS-1:0]   b_row;

    logic                            a_ready;
    logic                            a_valid;
    logic signed [IN_ROW_BITS-1:0]   a_row;

    logic                            c_valid;
    logic signed [OUT_ROW_BITS-1:0]  c_row;

    //-------------------------------------------------------------------------
    // DUT
    //-------------------------------------------------------------------------
    smm2_top_dbuf #(
        .TOP_N       (TOP_N),
        .SMM_R       (SMM_R),
        .BLOCK_N     (BLOCK_N),
        .TOP_G       (TOP_G),
        .DATA_W      (DATA_W),
        .TOP_MID_W   (TOP_MID_W),
        .SUB_IN_W    (SUB_IN_W),
        .SUB_MID_W   (SUB_MID_W),
        .ACC_W       (ACC_W),
        .SUB_OUT_W   (SUB_OUT_W),
        .OUT_W       (OUT_W)
    ) dut (
        .clk         (clk),
        .rst_n       (rst_n),

        .start       (start),
        .busy        (busy),
        .done        (done),
        .active_bank (active_bank),

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

    longint                   C_expected [TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_dut      [TOP_N][TOP_N];

    int c_idx;
    int errors;

    //-------------------------------------------------------------------------
    // Capture C rows from DUT
    //
    // Output packing:
    //   c_row = {C22_top, C21_top, C12_top, C11_top}
    // Each top quadrant is packed as {sub22, sub21, sub12, sub11}
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c_idx <= 0;
        end else if (c_valid && (c_idx < BLOCK_N)) begin

            for (int q = 0; q < 4; q++) begin : cap_q
                int tr, tc;

                tr = (q >= 2) ? 16 : 0;
                tc = (q == 1 || q == 3) ? 16 : 0;

                for (int s = 0; s < 4; s++) begin : cap_s
                    int sr, sc, row_base, col_base;

                    sr = (s >= 2) ? 1 : 0;
                    sc = (s == 1 || s == 3) ? 1 : 0;

                    row_base = tr + sr * BLOCK_N;
                    col_base = tc + sc * BLOCK_N;

                    for (int j = 0; j < BLOCK_N; j++) begin
                        C_dut[row_base + c_idx][col_base + j] <=
                            c_row[q*QUAD_OUT_BITS + s*SUB_OUT_BITS + j*OUT_W +: OUT_W];
                    end
                end
            end

            c_idx <= c_idx + 1;
        end
    end

    //-------------------------------------------------------------------------
    // Global watchdog
    //-------------------------------------------------------------------------
    initial begin
        #20_000_000;
        $display("TB ERROR: global watchdog timeout - simulation hung.");
        $finish;
    end

    //-------------------------------------------------------------------------
    // Reset
    //-------------------------------------------------------------------------
    task automatic reset_dut();
        begin
            rst_n    = 1'b0;
            start    = 1'b0;
            b_valid  = 1'b0;
            a_valid  = 1'b0;
            b_row    = '0;
            a_row    = '0;

            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++)
                    C_dut[r][c] = '0;

            repeat (10) @(posedge clk);
            rst_n = 1'b1;
            repeat (5) @(posedge clk);
        end
    endtask

    //-------------------------------------------------------------------------
    // Golden conventional model: C = A * B
    //-------------------------------------------------------------------------
    task automatic compute_expected();
        longint sum;
        begin
            for (int r = 0; r < TOP_N; r++) begin
                for (int c = 0; c < TOP_N; c++) begin
                    sum = 0;
                    for (int k = 0; k < TOP_N; k++) begin
                        sum = sum + (longint'(A[r][k]) * longint'(B[k][c]));
                    end
                    C_expected[r][c] = sum;
                end
            end
        end
    endtask

    //-------------------------------------------------------------------------
    // Test 1: A = identity, B deterministic
    //-------------------------------------------------------------------------
    task automatic init_identity_test();
        int val;
        begin
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++)
                    A[r][c] = (r == c) ? 16'sd1 : 16'sd0;

            for (int r = 0; r < TOP_N; r++) begin
                for (int c = 0; c < TOP_N; c++) begin
                    val      = (r * TOP_N + c) % 201;
                    B[r][c]  = val - 100;
                end
            end

            compute_expected();
        end
    endtask

    //-------------------------------------------------------------------------
    // Test 2: small deterministic ramp
    //-------------------------------------------------------------------------
    task automatic init_small_ramp_test();
        int val;
        begin
            for (int r = 0; r < TOP_N; r++) begin
                for (int c = 0; c < TOP_N; c++) begin
                    val     = ((r * TOP_N + c) % 17) - 8;
                    A[r][c] = val;

                    val     = ((r * TOP_N + c) % 13) - 6;
                    B[r][c] = val;
                end
            end

            compute_expected();
        end
    endtask

    //-------------------------------------------------------------------------
    // Packing helpers for A
    //-------------------------------------------------------------------------
    function automatic logic [SUB_IN_BITS-1:0] get_a_sub_row(
        input int q,   // top quadrant: 0=A11, 1=A12, 2=A21, 3=A22
        input int s,   // subquadrant: 0=sub11, 1=sub12, 2=sub21, 3=sub22
        input int i    // row index inside lowest 8x8 subblock
    );
        int tr, tc, sr, sc, row, col_base;
        begin
            tr = (q >= 2) ? 16 : 0;
            tc = (q == 1 || q == 3) ? 16 : 0;

            sr = (s >= 2) ? 1 : 0;
            sc = (s == 1 || s == 3) ? 1 : 0;

            row      = tr + sr * BLOCK_N + i;
            col_base = tc + sc * BLOCK_N;

            get_a_sub_row = '0;

            for (int j = 0; j < BLOCK_N; j++) begin
                get_a_sub_row[j*DATA_W +: DATA_W] = A[row][col_base + j];
            end
        end
    endfunction

    function automatic logic [QUAD_IN_BITS-1:0] get_a_quad_row(
        input int q,
        input int i
    );
        logic [SUB_IN_BITS-1:0] s11, s12, s21, s22;
        begin
            s11 = get_a_sub_row(q, 0, i);
            s12 = get_a_sub_row(q, 1, i);
            s21 = get_a_sub_row(q, 2, i);
            s22 = get_a_sub_row(q, 3, i);

            get_a_quad_row = {s22, s21, s12, s11};
        end
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_a_row(input int i);
        logic [QUAD_IN_BITS-1:0] q11, q12, q21, q22;
        begin
            q11 = get_a_quad_row(0, i);
            q12 = get_a_quad_row(1, i);
            q21 = get_a_quad_row(2, i);
            q22 = get_a_quad_row(3, i);

            get_a_row = {q22, q21, q12, q11};
        end
    endfunction

    //-------------------------------------------------------------------------
    // Packing helpers for B
    //-------------------------------------------------------------------------
    function automatic logic [SUB_IN_BITS-1:0] get_b_sub_row(
        input int q,
        input int s,
        input int i
    );
        int tr, tc, sr, sc, row, col_base;
        begin
            tr = (q >= 2) ? 16 : 0;
            tc = (q == 1 || q == 3) ? 16 : 0;

            sr = (s >= 2) ? 1 : 0;
            sc = (s == 1 || s == 3) ? 1 : 0;

            row      = tr + sr * BLOCK_N + i;
            col_base = tc + sc * BLOCK_N;

            get_b_sub_row = '0;

            for (int j = 0; j < BLOCK_N; j++) begin
                get_b_sub_row[j*DATA_W +: DATA_W] = B[row][col_base + j];
            end
        end
    endfunction

    function automatic logic [QUAD_IN_BITS-1:0] get_b_quad_row(
        input int q,
        input int i
    );
        logic [SUB_IN_BITS-1:0] s11, s12, s21, s22;
        begin
            s11 = get_b_sub_row(q, 0, i);
            s12 = get_b_sub_row(q, 1, i);
            s21 = get_b_sub_row(q, 2, i);
            s22 = get_b_sub_row(q, 3, i);

            get_b_quad_row = {s22, s21, s12, s11};
        end
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_b_row(input int i);
        logic [QUAD_IN_BITS-1:0] q11, q12, q21, q22;
        begin
            q11 = get_b_quad_row(0, i);
            q12 = get_b_quad_row(1, i);
            q21 = get_b_quad_row(2, i);
            q22 = get_b_quad_row(3, i);

            get_b_row = {q22, q21, q12, q11};
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
            while (!done && t < 200000) begin
                @(posedge clk);
                t = t + 1;
            end

            if (!done) begin
                $fatal(1, "TB ERROR: timeout waiting for done.");
            end

            repeat (10) @(posedge clk);
        end
    endtask

    //-------------------------------------------------------------------------
    // Compare results
    //-------------------------------------------------------------------------
    task automatic compare_results(input string test_name);
        longint dut_val;
        longint exp_val;
        begin
            errors = 0;

            if (c_idx != BLOCK_N) begin
                $display("TB ERROR: %s captured %0d C rows, expected %0d.",
                         test_name, c_idx, BLOCK_N);
                errors = errors + 1;
            end

            for (int r = 0; r < TOP_N; r++) begin
                for (int c = 0; c < TOP_N; c++) begin
                    dut_val = longint'(C_dut[r][c]);
                    exp_val = C_expected[r][c];

                    if (dut_val !== exp_val) begin
                        errors = errors + 1;

                        if (errors <= 20) begin
                            $display("TB ERROR: %s C[%0d][%0d] DUT=%0d EXP=%0d",
                                     test_name, r, c, dut_val, exp_val);
                        end
                    end
                end
            end

            $display("==================================================");
            if (errors == 0)
                $display("%s : PASS", test_name);
            else
                $display("%s : FAIL (%0d errors)", test_name, errors);
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

            @(negedge clk);
            start = 1'b1;
            @(negedge clk);
            start = 1'b0;

            $display("TB: [%s] start pulsed at %0t (active_bank=%0b)",
                     test_name, $time, active_bank);

            for (int i = 0; i < BLOCK_N; i++) begin
                send_b_row(get_b_row(i));
            end

            $display("TB: [%s] B tile sent at %0t (active_bank=%0b)",
                     test_name, $time, active_bank);

            for (int i = 0; i < BLOCK_N; i++) begin
                send_a_row(get_a_row(i));
            end

            $display("TB: [%s] A tile sent at %0t (active_bank=%0b)",
                     test_name, $time, active_bank);

            wait_done();

            $display("TB: [%s] done at %0t (active_bank=%0b)",
                     test_name, $time, active_bank);

            compare_results(test_name);
        end
    endtask

    //-------------------------------------------------------------------------
    // Main
    //-------------------------------------------------------------------------
    initial begin
        rst_n    = 1'b0;
        start    = 1'b0;
        b_valid  = 1'b0;
        a_valid  = 1'b0;
        b_row    = '0;
        a_row    = '0;

        init_identity_test();
        run_test("SMM2-DBUF Test 1: Identity A");

        init_small_ramp_test();
        run_test("SMM2-DBUF Test 2: Small ramp");

        $display("TB: simulation finished.");
        $finish;
    end

endmodule
