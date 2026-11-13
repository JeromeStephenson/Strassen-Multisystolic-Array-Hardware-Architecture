//=============================================================================
// tb_smm2_top_gemm.sv  (FIXED: full-rate drivers + steady-state metrics)
//=============================================================================
`timescale 1ns / 1ps

module tb_smm2_top_gemm;

    //-------------------------------------------------------------------------
    // Parameters
    //-------------------------------------------------------------------------
    localparam int TOP_N         = 32;
    localparam int SMM_R         = 2;
    localparam int BLOCK_N       = 8;
    localparam int TOP_G         = 4;
    localparam int NUM_TILES     = 16;   // FIX 2: enough tiles for steady state

    localparam int DATA_W        = 16;
    localparam int TOP_MID_W     = DATA_W + 1;
    localparam int SUB_IN_W      = TOP_MID_W;
    localparam int SUB_MID_W     = SUB_IN_W + 1;
    localparam int ACC_W         = 48;
    localparam int SUB_OUT_W     = ACC_W + 2;
    localparam int OUT_W         = SUB_OUT_W + 2;

    localparam int SUB_ELEMS     = BLOCK_N;
    localparam int QUAD_ELEMS    = 4 * BLOCK_N;
    localparam int SUB_IN_BITS   = SUB_ELEMS * DATA_W;
    localparam int QUAD_IN_BITS  = QUAD_ELEMS * DATA_W;
    localparam int IN_ROW_BITS   = 4 * QUAD_IN_BITS;
    localparam int SUB_OUT_BITS  = SUB_ELEMS * OUT_W;
    localparam int QUAD_OUT_BITS = QUAD_ELEMS * OUT_W;
    localparam int OUT_ROW_BITS  = 4 * QUAD_OUT_BITS;

    localparam int MXU_MULTS_ISO = 3136;              // 7^2 * 8^2  (Table I)
    localparam int MXU_MULTS_SYS = 3136 + 8*(4**2);   // +128 post-GEMM (Tables II/III) = 3264

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
        .TOP_N(TOP_N), .SMM_R(SMM_R), .BLOCK_N(BLOCK_N), .TOP_G(TOP_G),
        .DATA_W(DATA_W), .TOP_MID_W(TOP_MID_W), .SUB_IN_W(SUB_IN_W),
        .SUB_MID_W(SUB_MID_W), .ACC_W(ACC_W), .SUB_OUT_W(SUB_OUT_W), .OUT_W(OUT_W)
    ) dut (
        .clk(clk), .rst_n(rst_n), .start(start), .num_tiles(num_tiles),
        .busy(busy), .done(done), .active_bank(active_bank),
        .b_ready(b_ready), .b_valid(b_valid), .b_row(b_row),
        .a_ready(a_ready), .a_valid(a_valid), .a_row(a_row),
        .c_valid(c_valid), .c_row(c_row)
    );

    //-------------------------------------------------------------------------
    // Matrices, capture storage, metrics
    //-------------------------------------------------------------------------
    logic signed [DATA_W-1:0] A_mem      [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [DATA_W-1:0] B_mem      [0:NUM_TILES-1][TOP_N][TOP_N];
    longint                   C_expected [0:NUM_TILES-1][TOP_N][TOP_N];
    logic signed [OUT_W-1:0]  C_dut      [0:NUM_TILES-1][TOP_N][TOP_N];

    int c_count;
    int errors;

    longint start_cycle, end_cycle, total_cycles;
    longint total_mults;
    real    mults_per_cycle, mce, mse, throughput_per_dsp;

    // FIX 3: phase monitors
    int      a_seen, c_seen;
    longint  b_start_cyc [0:NUM_TILES-1];
    longint  b_end_cyc   [0:NUM_TILES-1];
    longint  a_start_cyc [0:NUM_TILES-1];
    longint  a_end_cyc   [0:NUM_TILES-1];
    longint  c_end_cyc   [0:NUM_TILES-1];

    //-------------------------------------------------------------------------
    // Capture C rows
    //-------------------------------------------------------------------------
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) c_count <= 0;
        else if (c_valid && (c_count < NUM_TILES*BLOCK_N)) begin
            automatic int t = c_count / BLOCK_N;
            automatic int i = c_count % BLOCK_N;
            for (int q = 0; q < 4; q++) begin
                int tr = (q >= 2) ? 16 : 0;
                int tc = (q == 1 || q == 3) ? 16 : 0;
                for (int s = 0; s < 4; s++) begin
                    int sr = (s >= 2) ? 1 : 0;
                    int sc = (s == 1 || s == 3) ? 1 : 0;
                    int row_base = tr + sr * BLOCK_N;
                    int col_base = tc + sc * BLOCK_N;
                    for (int j = 0; j < BLOCK_N; j++)
                        C_dut[t][row_base + i][col_base + j] <=
                            c_row[q*QUAD_OUT_BITS + s*SUB_OUT_BITS + j*OUT_W +: OUT_W];
                end
            end
            c_count <= c_count + 1;
        end
    end

    //-------------------------------------------------------------------------
    // FIX 3: handshake event monitor (timestamps of every tile phase)
    //-------------------------------------------------------------------------
    always @(posedge clk) begin
        if (a_valid && a_ready && (a_seen < NUM_TILES*BLOCK_N)) begin
            if (a_seen % BLOCK_N == 0)            a_start_cyc[a_seen/BLOCK_N] = $time/10;
            if (a_seen % BLOCK_N == BLOCK_N-1)    a_end_cyc  [a_seen/BLOCK_N] = $time/10;
            a_seen = a_seen + 1;
        end
        if (c_valid && (c_seen < NUM_TILES*BLOCK_N)) begin
            if (c_seen % BLOCK_N == BLOCK_N-1)    c_end_cyc  [c_seen/BLOCK_N] = $time/10;
            c_seen = c_seen + 1;
        end
    end

    initial begin #(NUM_TILES * 2_000_000); $display("TB ERROR: global watchdog timeout."); $finish; end

    //-------------------------------------------------------------------------
    // Packing helpers (unchanged)
    //-------------------------------------------------------------------------
    function automatic logic [SUB_IN_BITS-1:0] get_a_sub_row(input int t, input int q, input int s, input int i);
        int tr = (q >= 2) ? 16 : 0; int tc = (q == 1 || q == 3) ? 16 : 0;
        int sr = (s >= 2) ? 1 : 0;  int sc = (s == 1 || s == 3) ? 1 : 0;
        int row = tr + sr * BLOCK_N + i; int col_base = tc + sc * BLOCK_N;
        get_a_sub_row = '0;
        for (int j = 0; j < BLOCK_N; j++) get_a_sub_row[j*DATA_W +: DATA_W] = A_mem[t][row][col_base + j];
    endfunction

    function automatic logic [QUAD_IN_BITS-1:0] get_a_quad_row(input int t, input int q, input int i);
        get_a_quad_row = {get_a_sub_row(t,q,3,i), get_a_sub_row(t,q,2,i),
                          get_a_sub_row(t,q,1,i), get_a_sub_row(t,q,0,i)};
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_a_row(input int t, input int i);
        get_a_row = {get_a_quad_row(t,3,i), get_a_quad_row(t,2,i),
                     get_a_quad_row(t,1,i), get_a_quad_row(t,0,i)};
    endfunction

    function automatic logic [SUB_IN_BITS-1:0] get_b_sub_row(input int t, input int q, input int s, input int i);
        int tr = (q >= 2) ? 16 : 0; int tc = (q == 1 || q == 3) ? 16 : 0;
        int sr = (s >= 2) ? 1 : 0;  int sc = (s == 1 || s == 3) ? 1 : 0;
        int row = tr + sr * BLOCK_N + i; int col_base = tc + sc * BLOCK_N;
        get_b_sub_row = '0;
        for (int j = 0; j < BLOCK_N; j++) get_b_sub_row[j*DATA_W +: DATA_W] = B_mem[t][row][col_base + j];
    endfunction

    function automatic logic [QUAD_IN_BITS-1:0] get_b_quad_row(input int t, input int q, input int i);
        get_b_quad_row = {get_b_sub_row(t,q,3,i), get_b_sub_row(t,q,2,i),
                          get_b_sub_row(t,q,1,i), get_b_sub_row(t,q,0,i)};
    endfunction

    function automatic logic [IN_ROW_BITS-1:0] get_b_row(input int t, input int i);
        get_b_row = {get_b_quad_row(t,3,i), get_b_quad_row(t,2,i),
                     get_b_quad_row(t,1,i), get_b_quad_row(t,0,i)};
    endfunction

    //-------------------------------------------------------------------------
    // Stimulus / reference
    //-------------------------------------------------------------------------
    task automatic init_matrices();
        for (int t = 0; t < NUM_TILES; t++)
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    A_mem[t][r][c] = ((r * TOP_N + c + 7*t) % 17) - 8;
                    B_mem[t][r][c] = ((r * TOP_N + c + 5*t) % 13) - 6;
                end
    endtask

    task automatic compute_expected_all();
        longint sum;
        for (int t = 0; t < NUM_TILES; t++)
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    sum = 0;
                    for (int k = 0; k < TOP_N; k++)
                        sum = sum + (longint'(A_mem[t][r][k]) * longint'(B_mem[t][k][c]));
                    C_expected[t][r][c] = sum;
                end
    endtask

    task automatic compare_all();
        longint dut_val, exp_val;
        errors = 0;
        for (int t = 0; t < NUM_TILES; t++) begin
            int tile_err = 0;
            for (int r = 0; r < TOP_N; r++)
                for (int c = 0; c < TOP_N; c++) begin
                    dut_val = longint'(C_dut[t][r][c]); exp_val = C_expected[t][r][c];
                    if (dut_val !== exp_val) begin
                        tile_err++; errors++;
                        if (tile_err <= 5) $display("TB ERROR: tile%0d C[%0d][%0d] DUT=%0d EXP=%0d", t, r, c, dut_val, exp_val);
                    end
                end
            if (tile_err != 0) $display("SMM2-GEMM Tile %0d : FAIL (%0d errors)", t, tile_err);
        end
        if (errors == 0) $display("ALL %0d TILES PASS.", NUM_TILES);
    endtask

    //-------------------------------------------------------------------------
    // Main
    //-------------------------------------------------------------------------
    initial begin
        real ss_period, mpc_ss, mce_iso_ss, mce_sys_ss, mse_ss;
        real avg_bload, avg_aload, avg_drain_lat, avg_drain_wait;
        rst_n = 1'b0; start = 1'b0; num_tiles = NUM_TILES;
        b_valid = 1'b0; a_valid = 1'b0; b_row = '0; a_row = '0;
        a_seen = 0; c_seen = 0;
        init_matrices(); compute_expected_all();
        repeat (10) @(posedge clk); rst_n = 1'b1; repeat (5) @(posedge clk);
        @(negedge clk); start = 1'b1; start_cycle = $time / 10;
        @(negedge clk); start = 1'b0;

        begin
            automatic int tmo = 0;
            while (!done && tmo < 2_000_000) begin @(posedge clk); tmo++; end
            if (!done) $fatal(1, "TB ERROR: timeout waiting for done.");
            end_cycle = $time / 10;
        end
        repeat (10) @(posedge clk);
        compare_all();

        total_cycles = end_cycle - start_cycle;
        total_mults  = NUM_TILES * (TOP_N * TOP_N * TOP_N);

        //----- steady-state period: first-A-row of tile 0 -> first-A-row of last tile
        ss_period   = real'(a_start_cyc[NUM_TILES-1] - a_start_cyc[0]) / real'(NUM_TILES-1);
        mpc_ss      = real'(TOP_N*TOP_N*TOP_N) / ss_period;
        mce_iso_ss  = mpc_ss / real'(MXU_MULTS_ISO);
        mce_sys_ss  = mpc_ss / real'(MXU_MULTS_SYS);
        mse_ss      = mpc_ss / 1024.0;

        //----- phase averages over the steady region (tiles 1..N-2)
        avg_bload = 0; avg_aload = 0; avg_drain_lat = 0; avg_drain_wait = 0;
        for (int t = 1; t <= NUM_TILES-2; t++) begin
            avg_bload      += real'(b_end_cyc[t] - b_start_cyc[t]) + 1.0;
            avg_aload      += real'(a_end_cyc[t] - a_start_cyc[t]) + 1.0;
            avg_drain_lat  += real'(c_end_cyc[t] - a_end_cyc[t])   + 1.0;
            avg_drain_wait += real'(a_start_cyc[t+1] - c_end_cyc[t]) - 1.0;
        end
        avg_bload      /= real'(NUM_TILES-2);
        avg_aload      /= real'(NUM_TILES-2);
        avg_drain_lat  /= real'(NUM_TILES-2);
        avg_drain_wait /= real'(NUM_TILES-2);

        //----- whole-run (legacy) metrics
        mults_per_cycle    = real'(total_mults) / real'(total_cycles);
        mce                = mults_per_cycle / real'(MXU_MULTS_ISO);
        mse                = mults_per_cycle / 1024.0;
        throughput_per_dsp = mults_per_cycle / real'(MXU_MULTS_ISO);

        $display("==================================================");
        $display(" PAPER METRICS (SMM2, r=2, n=32, %0d tiles)", NUM_TILES);
        $display("==================================================");
        $display(" Total cycles (run)         : %0d", total_cycles);
        $display(" STEADY-STATE cycles/tile   : %0.2f  (roof = 8.00)", ss_period);
        $display(" --- phase breakdown (cycles) ---");
        $display("   B load / tile            : %0.2f  (hidden if < period)", avg_bload);
        $display("   A stream / tile          : %0.2f  (roof = 8)", avg_aload);
        $display("   drain latency            : %0.2f  (pipeline depth)", avg_drain_lat);
        $display("   IDLE gap after drain     : %0.2f  (roof <= 0 => overlapped)", avg_drain_wait);
        $display(" --- steady-state metrics ---");
        $display(" Mults / clock cycle        : %0.2f  (roof = 4096)", mpc_ss);
        $display(" MCE isolated  (/3136)      : %0.4f  (roof = 1.306)", mce_iso_ss);
        $display(" MCE system    (/3264)      : %0.4f  (Tables II/III)", mce_sys_ss);
        $display(" MSE           (/1024)      : %0.4f  (roof = 4.0)", mse_ss);
        $display(" --- whole-run metrics (incl. prologue) ---");
        $display(" Mults / clock cycle        : %0.2f", mults_per_cycle);
        $display(" MCE (/3136)                : %0.4f", mce);
        $display(" Throughput / DSP           : %0.4f mults/cycle/mult", throughput_per_dsp);
        $display(" DSP saving vs MM2          : %0.1f%%", ((4096.0/3136.0)-1.0)*100.0);
        $display("==================================================");
        $display("TB: simulation finished (errors=%0d).", errors);
        $finish;
    end

    //-------------------------------------------------------------------------
    // FIX 1: FULL-RATE B driver (1 row per cycle, valid held until handshake)
    //-------------------------------------------------------------------------
    initial begin
        b_valid = 1'b0; b_row = '0;
        wait (rst_n === 1'b1);
        for (int t = 0; t < NUM_TILES; t++) begin
            for (int i = 0; i < BLOCK_N; i++) begin
                @(negedge clk);
                b_valid = 1'b1;
                b_row   = get_b_row(t, i);
                @(posedge clk);
                while (!b_ready) @(posedge clk);   // hold if stalled
                if (i == 0) b_start_cyc[t] = $time/10 - 1;
                if (i == BLOCK_N-1) b_end_cyc[t] = $time/10 - 1;
            end
        end
        @(negedge clk); b_valid = 1'b0; b_row = '0;
    end

    //-------------------------------------------------------------------------
    // FIX 1: FULL-RATE A driver
    //-------------------------------------------------------------------------
    initial begin
        a_valid = 1'b0; a_row = '0;
        wait (rst_n === 1'b1);
        for (int t = 0; t < NUM_TILES; t++) begin
            for (int i = 0; i < BLOCK_N; i++) begin
                @(negedge clk);
                a_valid = 1'b1;
                a_row   = get_a_row(t, i);
                @(posedge clk);
                while (!a_ready) @(posedge clk);   // hold if stalled
            end
        end
        @(negedge clk); a_valid = 1'b0; a_row = '0;
    end

endmodule
