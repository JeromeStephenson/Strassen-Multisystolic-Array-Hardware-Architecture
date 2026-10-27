//=============================================================================
// smm1_top.sv
//
// Top-level SMM1 datapath and controller for one 16x16 matrix tile.
//
// Operation:
//   1. Assert start.
//   2. Feed 8 B block-row vectors while b_ready is high.
//   3. Feed 8 A block-row vectors while a_ready is high.
//   4. Receive 8 C block-row vectors on c_valid/c_row.
//
// Input packing:
//   a_row = {A11 row, A12 row, A21 row, A22 row}
//   b_row = {B11 row, B12 row, B21 row, B22 row}
//
// Output packing:
//   c_row = {C11 row, C12 row, C21 row, C22 row}
//=============================================================================

(* keep_hierarchy = "yes" *)
module smm1_top #(
    parameter int TOP_N    = 16,
    parameter int SMM_R    = 1,
    parameter int BLOCK_N  = TOP_N / (1 << SMM_R), // 8 for r=1
    parameter int DATA_W   = 16,
    parameter int MID_W    = DATA_W + SMM_R,       // 17
    parameter int ACC_W    = 48,
    parameter int OUT_W    = ACC_W + 2             // 50
) (
    input  logic                                    clk,
    input  logic                                    rst_n,

    input  logic                                    start,
    output logic                                    busy,
    output logic                                    done,

    output logic                                    b_ready,
    input  logic                                    b_valid,
    input  logic signed [4*BLOCK_N*DATA_W-1:0]      b_row,

    output logic                                    a_ready,
    input  logic                                    a_valid,
    input  logic signed [4*BLOCK_N*DATA_W-1:0]      a_row,

    output logic                                    c_valid,
    output logic signed [4*BLOCK_N*OUT_W-1:0]       c_row
);

    //-------------------------------------------------------------------------
    // Local parameters
    //-------------------------------------------------------------------------
    localparam int CNT_W     = $clog2(BLOCK_N) + 1;
    localparam int ROW_IDX_W = (BLOCK_N > 1) ? $clog2(BLOCK_N) : 1;

    //-------------------------------------------------------------------------
    // FSM
    //-------------------------------------------------------------------------
    typedef enum logic [1:0] {
        ST_IDLE,
        ST_LOAD_B,
        ST_COMPUTE,
        ST_DONE
    } state_e;

    state_e state, state_next;

    //-------------------------------------------------------------------------
    // Counters
    //-------------------------------------------------------------------------
    logic [CNT_W-1:0] b_cnt;
    logic [CNT_W-1:0] s_load_cnt;
    logic [CNT_W-1:0] a_cnt;
    logic [CNT_W-1:0] c_cnt;

    //-------------------------------------------------------------------------
    // Delayed load row (aligns with s_valid, one cycle after pre_b_valid)
    //-------------------------------------------------------------------------
    logic [ROW_IDX_W-1:0] load_row_d;

    //-------------------------------------------------------------------------
    // Internal signals
    //-------------------------------------------------------------------------
    logic pre_a_valid;
    logic pre_b_valid;
    logic t_valid;
    logic s_valid;
    logic load_en;

    logic signed [BLOCK_N*MID_W-1:0] t1, t2, t3, t4, t5, t6, t7;
    logic signed [BLOCK_N*MID_W-1:0] s1, s2, s3, s4, s5, s6, s7;

    logic q_valid;
    logic signed [BLOCK_N*ACC_W-1:0] q1, q2, q3, q4, q5, q6, q7;

    logic c_valid_post;
    logic signed [BLOCK_N*OUT_W-1:0] c11_row, c12_row, c21_row, c22_row;

    //-------------------------------------------------------------------------
    // Input handshakes
    //-------------------------------------------------------------------------
    assign b_ready = (state == ST_LOAD_B) && (b_cnt < CNT_W'(BLOCK_N));
    assign a_ready = (state == ST_COMPUTE) && (a_cnt < CNT_W'(BLOCK_N));

    assign pre_b_valid = b_ready & b_valid;
    assign pre_a_valid = a_ready & a_valid;

    //-------------------------------------------------------------------------
    // Load row delay
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            load_row_d <= '0;
        end else begin
            if (pre_b_valid) begin
                load_row_d <= b_cnt[ROW_IDX_W-1:0];
            end
        end
    end

    //-------------------------------------------------------------------------
    // Load enable (s_valid is pre_b_valid delayed by one cycle)
    //-------------------------------------------------------------------------
    assign load_en = (state == ST_LOAD_B) & s_valid;

    //-------------------------------------------------------------------------
    // Counters
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            b_cnt      <= '0;
            s_load_cnt <= '0;
            a_cnt      <= '0;
            c_cnt      <= '0;
        end else begin
            case (state)
                ST_IDLE: begin
                    if (start) begin
                        b_cnt      <= '0;
                        s_load_cnt <= '0;
                        a_cnt      <= '0;
                        c_cnt      <= '0;
                    end
                end
                ST_LOAD_B: begin
                    if (pre_b_valid) b_cnt <= b_cnt + CNT_W'(1);
                    if (load_en)     s_load_cnt <= s_load_cnt + CNT_W'(1);
                end
                ST_COMPUTE: begin
                    if (pre_a_valid) a_cnt <= a_cnt + CNT_W'(1);
                    if (c_valid)     c_cnt <= c_cnt + CNT_W'(1);
                end
                default: ;
            endcase
        end
    end

    //-------------------------------------------------------------------------
    // FSM next-state logic
    //-------------------------------------------------------------------------
    always_comb begin
        state_next = state;
        case (state)
            ST_IDLE:    if (start) state_next = ST_LOAD_B;
            ST_LOAD_B:  if (s_load_cnt == CNT_W'(BLOCK_N)) state_next = ST_COMPUTE;
            ST_COMPUTE: if (c_cnt == CNT_W'(BLOCK_N)) state_next = ST_DONE;
            ST_DONE:    state_next = ST_IDLE;
            default:    state_next = ST_IDLE;
        endcase
    end

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) state <= ST_IDLE;
        else        state <= state_next;
    end

    assign busy = (state == ST_LOAD_B) || (state == ST_COMPUTE);
    assign done = (state == ST_DONE);

    //-------------------------------------------------------------------------
    // Strassen pre-addition network
    //-------------------------------------------------------------------------
    strassen_pre #(
        .N(BLOCK_N), .W(DATA_W), .MW(MID_W)
    ) u_strassen_pre (
        .clk(clk), .rst_n(rst_n),
        .a_valid(pre_a_valid), .b_valid(pre_b_valid),
        .a_row(a_row), .b_row(b_row),
        .t1(t1), .t2(t2), .t3(t3), .t4(t4), .t5(t5), .t6(t6), .t7(t7),
        .s1(s1), .s2(s2), .s3(s3), .s4(s4), .s5(s5), .s6(s6), .s7(s7),
        .t_valid(t_valid), .s_valid(s_valid)
    );

    //-------------------------------------------------------------------------
    // Seven Strassen MM cores
    //-------------------------------------------------------------------------
    smm1_cores #(
        .N(BLOCK_N), .MW(MID_W), .ACC_W(ACC_W)
    ) u_smm1_cores (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row_d),
        .s1(s1), .s2(s2), .s3(s3), .s4(s4), .s5(s5), .s6(s6), .s7(s7),
        .t_valid(t_valid),
        .t1(t1), .t2(t2), .t3(t3), .t4(t4), .t5(t5), .t6(t6), .t7(t7),
        .q_valid(q_valid),
        .q1(q1), .q2(q2), .q3(q3), .q4(q4), .q5(q5), .q6(q6), .q7(q7)
    );

    //-------------------------------------------------------------------------
    // Strassen post-addition network
    //-------------------------------------------------------------------------
    strassen_post #(
        .N(BLOCK_N), .QW(ACC_W), .CW(OUT_W)
    ) u_strassen_post (
        .clk(clk), .rst_n(rst_n),
        .in_valid(q_valid),
        .q1(q1), .q2(q2), .q3(q3), .q4(q4), .q5(q5), .q6(q6), .q7(q7),
        .c11(c11_row), .c12(c12_row), .c21(c21_row), .c22(c22_row),
        .out_valid(c_valid_post)
    );

    //-------------------------------------------------------------------------
    // Output C stream
    //-------------------------------------------------------------------------
    assign c_valid = c_valid_post & (state == ST_COMPUTE);
    assign c_row   = {c22_row, c21_row, c12_row, c11_row};

endmodule
