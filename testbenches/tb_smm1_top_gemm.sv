//=============================================================================
// tb_smm1_top_gemm.sv
//
// Testbench for smm1_top_gemm.sv (paper-style overlapped multi-tile top).
//
// num_tiles = 2:
//   C_tile_0 = A_tile_0 * B_tile_0
//   C_tile_1 = A_tile_1 * B_tile_1
//
// B driver and A driver run concurrently, so B_tile_1 is loaded into the
// shadow buffer WHILE A_tile_0 is being multiplied (paper-style overlap).
//=============================================================================

`timescale 1ns / 1ps

module tb_smm1_top_gemm;

    //-------------------------------------------------------------------------
    // Parameters
    //-------------------------------------------------------------------------
    localparam int TOP_N     = 16;
    localparam int SMM_R     = 1;
    localparam int BLOCK_N   = 8;
    localparam int DATA_W    = 16;
    localparam int MID_W     = DATA_W + SMM_R;
    localparam int ACC_W     = 48;
    localparam int OUT_W     = ACC_W + 2;
    localparam int NUM_TILES = 2;

    localparam int IN_ROW_BITS   = 4 * BLOCK_N * DATA_W;
    localparam int OUT_ROW_BITS  = 4 * BLOCK_N * OUT_W;
    localparam int QUAD_IN_BITS  = BLOCK_N * DATA_W;
    localparam int QUAD_OUT_BITS = BLOCK_N * OUT_W;

    //-------------------------------------------------------------------------
    // Clock
    //-------------------------------------------------------------------------
    logic clk = 1'b0;
    always #5 clk = ~clk;

    //-------------------------------------------------------------------------
    // DUT interface
    //-------------------------------------------------------------------------
    logic                       rst_n;
    logic                       start;
    logic [15:0]                num_tiles;
    logic                       busy;
    logic                       done;
    logic                       active_bank;
    logic                       b_ready;
    logic                       b_valid;
    logic [IN_ROW_BITS-1:0]     b_row;
    logic                       a_ready;
    logic                       a_valid;
    logic [IN_ROW_BITS-1:0]     a_row;
    logic                       c_valid;
    logic [OUT_ROW_BITS-1:0]    c_row;

    smm1_top_gemm #(
        .TOP_N   (TOP_N),
        .SMM_R   (SMM_R),
        .BLOCK_N (BLOCK_N),
        .DATA_W  (DATA_W),
        .MID_W   (MID_W),
        .ACC_W   (ACC_W),
        .OUT_W   (OUT_W)
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
    logic signed [DATA_W-1:0] A_mem [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [DATA_W-1:0] B_mem [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_expected [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_dut      [0:NUM_TILES-1][TOP_N][TOP_N];

    int c_count;
    int errors;

    //-------------------------------------------------------------------------
    // Capture C rows (16 total: 8 for tile 0, 8 for tile 1)
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c_count <= 0;
        end else if (c_valid && (c_count < NUM_TILES*BLOCK_N)) begin
            for (int j = 0; j < BLOCK_N; j++) begin
                C_dut[c_count/BLOCK_N][c_count%BLOCK_N][j] <=
                    c_row[0*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
                C_dut[c_count/BLOCK_N][c_count%BLOCK_N][BLOCK_N+j] <=
                    c_row[1*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
                C_dut[c_count/BLOCK_N][BLOCK_N+(c_count%BLOCK_N)][j] <=
                    c_row[2*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
                C_dut[c_count/BLOCK_N][BLOCK_N+(c_count%BLOCK_N)][BLOCK_N+j] <=
                    c_row[3*QUAD_OUT_BITS + j*OUT_W +: OUT_W];
            end
            c_count <= c_count + 1;
        end
    end

    //-------------------------------------------------------------------------
    // Watchdog
    //-------------------------------------------------------------------------
    initial begin
        #4_000_000;
        $display("TB ERROR: global watchdog timeout.");
        $finish;
    end

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
    // Packing functions
    //-------------------------------------------------------------------------
    function automatic logic [IN_ROW_BITS-1:0] get_a_row(input int t, input int i);
        logic [QUAD_IN_BITS-1:0] a11, a12, a21, a22;
        begin
            a11 = '0; a12 = '0; a21 = '0; a22 = '0;
            for (int j = 0; j < BLOCK_N; j++) begin
                a11[j*DATA_W +: DATA_W] = A_mem[t][i][j];
                a12[j*DATA_W +: DATA_W] = A_mem[t][i][BLOCK_N+j];
                a21[j*DATA_W +: DATA_W] = A_mem[t][BLOCK_N+i][j];
                a22[j*DATA_W +: DATA_W] = A_mem[t][BLOCK_N+i][BLOCK_N+j];
            end
            get_a_row = {a22, a21, a12, a11};
        end
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_b_row(input int t, input int i);
        logic [QUAD_IN_BITS-1:0] b11, b12, b21, b22;
        begin
            b11 = '0; b12 = '0; b21 = '0; b22 = '0;
            for (int j = 0; j < BLOCK_N; j++) begin
                b11[j*DATA_W +: DATA_W] = B_mem[t][i][j];
                b12[j*DATA_W +: DATA_W] = B_mem[t][i][BLOCK_N+j];
                b21[j*DATA_W +: DATA_W] = B_mem[t][BLOCK_N+i][j];
                b22[j*DATA_W +: DATA_W] = B_mem[t][BLOCK_N+i][BLOCK_N+j];
            end
            get_b_row = {b22, b21, b12, b11};
        end
    endfunction

    //-------------------------------------------------------------------------
    // Initialize matrices (deterministic, different per tile)
    //-------------------------------------------------------------------------
    task automatic init_matrices();
        int va; int vb;
        begin
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    // Tile 0
                    va = ((r*TOP_N + c) % 17) - 8;
                    vb = ((r*TOP_N + c) % 13) - 6;
                    A_mem[0][r][c] = va;
                    B_mem[0][r][c] = vb;
                    // Tile 1 (different patterns, includes negatives)
                    va = ((r*TOP_N + c) % 19) - 9;
                    vb = ((r*TOP_N + c) % 23) - 11;
                    A_mem[1][r][c] = va;
                    B_mem[1][r][c] = vb;
                end
        end
    endtask

    //-------------------------------------------------------------------------
    // Golden model
    //-------------------------------------------------------------------------
    task automatic compute_expected_all();
        longint sum;
        begin
            for (int t = 0; t < NUM_TILES; t++)
                for (int r = 0; r < TOP_N; r++)
                    for (int c = 0; c < TOP_N; c++) begin
                        sum = 0;
                        for (int k = 0; k < TOP_N; k++)
                            sum = sum + (longint'(A_mem[t][r][k]) *
                                         longint'(B_mem[t][k][c]));
                        C_expected[t][r][c] = sum;
                    end
        end
    endtask

    //-------------------------------------------------------------------------
    // Compare
    //-------------------------------------------------------------------------
    task automatic compare_all();
        begin
            errors = 0;
            for (int t = 0; t < NUM_TILES; t++) begin
                int tile_err = 0;
                for (int r = 0; r < TOP_N; r++)
                    for (int c = 0; c < TOP_N; c++)
                        if (C_dut[t][r][c] !== C_expected[t][r][c]) begin
                            tile_err = tile_err + 1;
                            if (tile_err <= 10)
                                $display("TB ERROR: tile%0d C[%0d][%0d] DUT=%0d EXP=%0d",
                                         t, r, c, C_dut[t][r][c], C_expected[t][r][c]);
                        end
                if (tile_err == 0)
                    $display("==================================================");
                    if (tile_err == 0) $display("Tile %0d : PASS", t);
                    else               $display("Tile %0d : FAIL (%0d errors)", t, tile_err);
                    if (tile_err == 0)
                    $display("==================================================");
                errors = errors + tile_err;
            end
        end
    endtask

    //-------------------------------------------------------------------------
    // B driver (concurrent): streams B0 then B1.
    // B1 can only be accepted during tile-0 compute -> demonstrates overlap.
    //-------------------------------------------------------------------------
    initial begin
        b_valid = 1'b0; b_row = '0;
        wait (rst_n === 1'b1);
        for (int t = 0; t < NUM_TILES; t++) begin
            for (int i = 0; i < BLOCK_N; i++)
                send_b_row(get_b_row(t, i));
            $display("TB: B tile %0d fully sent at %0t", t, $time);
        end
    end

    //-------------------------------------------------------------------------
    // A driver (concurrent): streams A0 then A1.
    //-------------------------------------------------------------------------
    initial begin
        a_valid = 1'b0; a_row = '0;
        wait (rst_n === 1'b1);
        for (int t = 0; t < NUM_TILES; t++) begin
            for (int i = 0; i < BLOCK_N; i++)
                send_a_row(get_a_row(t, i));
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
            while (!done && tmo < 40000) begin
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
        $display("TB: simulation finished.");
        $finish;
    end

endmodule





