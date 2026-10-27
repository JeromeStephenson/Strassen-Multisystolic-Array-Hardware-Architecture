//=============================================================================
// tb_smm2_top_gemm.sv
//
// Testbench for smm2_top_gemm.sv (paper-style overlapped multi-tile SMM2).
//
//  num_tiles = 2:
//   C_tile_0 = A_tile_0 * B_tile_0
//   C_tile_1 = A_tile_1 * B_tile_1
//
// The B driver and A driver run concurrently, so B_tile_1 is loaded into
// the shadow buffer WHILE A_tile_0 is being multiplied (paper-style overlap).
//
// Uses the same recursive 32x32 packing as tb_smm2_top_dbuf.sv.
//=============================================================================

`timescale 1ns / 1ps

module tb_smm2_top_gemm;

    //-------------------------------------------------------------------------
    // Parameters (must match DUT)
    //-------------------------------------------------------------------------
    localparam int TOP_N         = 32;
    localparam int SMM_R         = 2;
    localparam int BLOCK_N       = 8;
    localparam int TOP_G         = 4;
    localparam int NUM_TILES     = 2;

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
    logic [15:0]                     num_tiles;
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
    smm2_top_gemm #(
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
        .num_tiles   (num_tiles),
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
    // Per-tile matrices
    //-------------------------------------------------------------------------
    logic signed [DATA_W-1:0] A_mem      [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [DATA_W-1:0] B_mem      [0:NUM_TILES-1][TOP_N][TOP_N];
    longint                   C_expected [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_dut      [0:NUM_TILES-1][TOP_N][TOP_N];

    int c_count;
    int errors;

    //-------------------------------------------------------------------------
    // Capture C rows from DUT (recursive unpack)
    //
    // c_count ranges 0 .. NUM_TILES*BLOCK_N-1
    //   tile            = c_count / BLOCK_N
    //   subblock row    = c_count % BLOCK_N
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c_count <= 0;
        end else if (c_valid && (c_count < NUM_TILES*BLOCK_N)) begin
            automatic int t = c_count / BLOCK_N;
            automatic int i = c_count % BLOCK_N;

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
                        C_dut[t][row_base + i][col_base + j] <=
                            c_row[q*QUAD_OUT_BITS + s*SUB_OUT_BITS + j*OUT_W +: OUT_W];
                    end
                end
            end

            c_count <= c_count + 1;
        end
    end

    //-------------------------------------------------------------------------
    // Global watchdog
    //-------------------------------------------------------------------------
    initial begin
        #40_000_000;
        $display("TB ERROR: global watchdog timeout - simulation hung.");
        $finish;
    end

    //-------------------------------------------------------------------------
    // Packing helpers for A (tile t, row i)
    //-------------------------------------------------------------------------
    function automatic logic [SUB_IN_BITS-1:0] get_a_sub_row(
        input int t,
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
                get_a_sub_row[j*DATA_W +: DATA_W] = A_mem[t][row][col_base + j];
            end
        end
    endfunction

    function automatic logic [QUAD_IN_BITS-1:0] get_a_quad_row(
        input int t,
        input int q,
        input int i
    );
        logic [SUB_IN_BITS-1:0] s11, s12, s21, s22;
        begin
            s11 = get_a_sub_row(t, q, 0, i);
            s12 = get_a_sub_row(t, q, 1, i);
            s21 = get_a_sub_row(t, q, 2, i);
            s22 = get_a_sub_row(t, q, 3, i);

            get_a_quad_row = {s22, s21, s12, s11};
        end
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_a_row(
        input int t,
        input int i
    );
        logic [QUAD_IN_BITS-1:0] q11, q12, q21, q22;
        begin
            q11 = get_a_quad_row(t, 0, i);
            q12 = get_a_quad_row(t, 1, i);
            q21 = get_a_quad_row(t, 2, i);
            q22 = get_a_quad_row(t, 3, i);

            get_a_row = {q22, q21, q12, q11};
        end
    endfunction

    //-------------------------------------------------------------------------
    // Packing helpers for B (tile t, row i)
    //-------------------------------------------------------------------------
    function automatic logic [SUB_IN_BITS-1:0] get_b_sub_row(
        input int t,
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
                get_b_sub_row[j*DATA_W +: DATA_W] = B_mem[t][row][col_base + j];
            end
        end
    endfunction

    function automatic logic [QUAD_IN_BITS-1:0] get_b_quad_row(
        input int t,
        input int q,
        input int i
    );
        logic [SUB_IN_BITS-1:0] s11, s12, s21, s22;
        begin
            s11 = get_b_sub_row(t, q, 0, i);
            s12 = get_b_sub_row(t, q, 1, i);
            s21 = get_b_sub_row(t, q, 2, i);
            s22 = get_b_sub_row(t, q, 3, i);

            get_b_quad_row = {s22, s21, s12, s11};
        end
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_b_row(
        input int t,
        input int i
    );
        logic [QUAD_IN_BITS-1:0] q11, q12, q21, q22;
        begin
            q11 = get_b_quad_row(t, 0, i);
            q12 = get_b_quad_row(t, 1, i);
            q21 = get_b_quad_row(t, 2, i);
            q22 = get_b_quad_row(t, 3, i);

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
    // Initialize matrices (deterministic, different per tile)
    //-------------------------------------------------------------------------
    task automatic init_matrices();
        int va;
        int vb;
        begin
            for (int r = 0; r < TOP_N; r++) begin
                for (int c = 0; c < TOP_N; c++) begin
                    // Tile 0
                    va = ((r * TOP_N + c) % 17) - 8;
                    vb = ((r * TOP_N + c) % 13) - 6;
                    A_mem[0][r][c] = va;
                    B_mem[0][r][c] = vb;

                    // Tile 1 (different patterns, includes negatives)
                    va = ((r * TOP_N + c) % 19) - 9;
                    vb = ((r * TOP_N + c) % 23) - 11;
                    A_mem[1][r][c] = va;
                    B_mem[1][r][c] = vb;
                end
            end
        end
    endtask

    //-------------------------------------------------------------------------
    // Golden conventional model: C = A * B for each tile
    //-------------------------------------------------------------------------
    task automatic compute_expected_all();
        longint sum;
        begin
            for (int t = 0; t < NUM_TILES; t++) begin
                for (int r = 0; r < TOP_N; r++) begin
                    for (int c = 0; c < TOP_N; c++) begin
                        sum = 0;
                        for (int k = 0; k < TOP_N; k++) begin
                            sum = sum + (longint'(A_mem[t][r][k]) * longint'(B_mem[t][k][c]));
                        end
                        C_expected[t][r][c] = sum;
                    end
                end
            end
        end
    endtask

    //-------------------------------------------------------------------------
    // Compare
    //-------------------------------------------------------------------------
    task automatic compare_all();
        longint dut_val;
        longint exp_val;
        begin
            errors = 0;

            for (int t = 0; t < NUM_TILES; t++) begin
                int tile_err;
                tile_err = 0;

                for (int r = 0; r < TOP_N; r++) begin
                    for (int c = 0; c < TOP_N; c++) begin
                        dut_val = longint'(C_dut[t][r][c]);
                        exp_val = C_expected[t][r][c];

                        if (dut_val !== exp_val) begin
                            tile_err = tile_err + 1;
                            errors   = errors + 1;

                            if (tile_err <= 10) begin
                                $display("TB ERROR: tile%0d C[%0d][%0d] DUT=%0d EXP=%0d",
                                         t, r, c, dut_val, exp_val);
                            end
                        end
                    end
                end

                $display("==================================================");
                if (tile_err == 0)
                    $display("SMM2-GEMM Tile %0d : PASS", t);
                else
                    $display("SMM2-GEMM Tile %0d : FAIL (%0d errors)", t, tile_err);
                $display("==================================================");
            end
        end
    endtask

    //-------------------------------------------------------------------------
    // B driver (concurrent): streams B0 then B1.
    // B1 can only be accepted during tile-0 compute -> demonstrates overlap.
    //-------------------------------------------------------------------------
    initial begin
        b_valid = 1'b0;
        b_row   = '0;

        wait (rst_n === 1'b1);

        for (int t = 0; t < NUM_TILES; t++) begin
            for (int i = 0; i < BLOCK_N; i++) begin
                send_b_row(get_b_row(t, i));
            end
            $display("TB: B tile %0d fully sent at %0t", t, $time);
        end
    end

    //-------------------------------------------------------------------------
    // A driver (concurrent): streams A0 then A1.
    //-------------------------------------------------------------------------
    initial begin
        a_valid = 1'b0;
        a_row   = '0;

        wait (rst_n === 1'b1);

        for (int t = 0; t < NUM_TILES; t++) begin
            for (int i = 0; i < BLOCK_N; i++) begin
                send_a_row(get_a_row(t, i));
            end
            $display("TB: A tile %0d fully sent at %0t", t, $time);
        end
    end

    //-------------------------------------------------------------------------
    // Main
    //-------------------------------------------------------------------------
    initial begin
        rst_n     = 1'b0;
        start     = 1'b0;
        num_tiles = NUM_TILES;

        init_matrices();
        compute_expected_all();

        repeat (10) @(posedge clk);
        rst_n = 1'b1;
        repeat (5) @(posedge clk);

        @(negedge clk); start = 1'b1;
        @(negedge clk); start = 1'b0;

        $display("TB: start pulsed at %0t (num_tiles=%0d)", $time, num_tiles);

        begin
            automatic int tmo = 0;
            while (!done && tmo < 400000) begin
                @(posedge clk);
                tmo = tmo + 1;
            end
            if (!done) $fatal(1, "TB ERROR: timeout waiting for done.");
        end

        $display("TB: done at %0t", $time);

        repeat (10) @(posedge clk);

        compare_all();

        if (errors == 0)
            $display("ALL TILES PASS - paper-style overlap verified.");
        else
            $display("TOTAL ERRORS: %0d", errors);

        $display("TB: simulation finished.");
        $finish;
    end

endmodule
