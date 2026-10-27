//=============================================================================
// mm1_top.sv
//
// Baseline conventional multisystolic top (MM1) for one 16x16 tile.
// Used for the DSP comparison against smm1_top / smm1_top_dbuf.
//=============================================================================

(* keep_hierarchy = "yes" *)
module mm1_top #(
    parameter int TOP_N    = 16,
    parameter int SMM_R    = 1,
    parameter int BLOCK_N  = TOP_N / (1 << SMM_R),
    parameter int DATA_W   = 16,
    parameter int ACC_W    = 48,
    parameter int OUT_W    = ACC_W + 2
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

    localparam int CNT_W     = $clog2(BLOCK_N) + 1;
    localparam int ROW_IDX_W = (BLOCK_N > 1) ? $clog2(BLOCK_N) : 1;
    localparam int QW        = BLOCK_N * DATA_W;
    localparam int CW        = OUT_W;

    typedef enum logic [1:0] {ST_IDLE, ST_LOAD_B, ST_COMPUTE, ST_DONE} state_e;
    state_e state, state_next;

    logic [CNT_W-1:0] b_cnt, a_cnt, c_cnt;

    // Unpacked quadrants (no pre-addition in baseline)
    logic signed [QW-1:0] a11, a12, a21, a22;
    logic signed [QW-1:0] b11, b12, b21, b22;

    assign a11 = a_row[0*QW +: QW];
    assign a12 = a_row[1*QW +: QW];
    assign a21 = a_row[2*QW +: QW];
    assign a22 = a_row[3*QW +: QW];

    assign b11 = b_row[0*QW +: QW];
    assign b12 = b_row[1*QW +: QW];
    assign b21 = b_row[2*QW +: QW];
    assign b22 = b_row[3*QW +: QW];

    logic pre_a_valid, pre_b_valid, t_valid, q_valid;
    logic signed [BLOCK_N*ACC_W-1:0] q1, q2, q3, q4, q5, q6, q7, q8;

    //-------------------------------------------------------------------------
    // Handshakes (no pre-adder latency in baseline)
    //-------------------------------------------------------------------------
    assign b_ready = (state == ST_LOAD_B) && (b_cnt < CNT_W'(BLOCK_N));
    assign a_ready = (state == ST_COMPUTE) && (a_cnt < CNT_W'(BLOCK_N));
    assign pre_b_valid = b_ready & b_valid;
    assign pre_a_valid = a_ready & a_valid;
    assign t_valid     = pre_a_valid;

    //-------------------------------------------------------------------------
    // Counters / FSM
    //-------------------------------------------------------------------------
    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            b_cnt <= '0; a_cnt <= '0; c_cnt <= '0;
        end else begin
            case (state)
                ST_IDLE:    if (start) begin b_cnt <= '0; a_cnt <= '0; c_cnt <= '0; end
                ST_LOAD_B:  if (pre_b_valid) b_cnt <= b_cnt + CNT_W'(1);
                ST_COMPUTE: begin
                    if (pre_a_valid) a_cnt <= a_cnt + CNT_W'(1);
                    if (c_valid)     c_cnt <= c_cnt + CNT_W'(1);
                end
                default: ;
            endcase
        end
    end

    always_comb begin
        state_next = state;
        case (state)
            ST_IDLE:    if (start) state_next = ST_LOAD_B;
            ST_LOAD_B:  if (b_cnt == CNT_W'(BLOCK_N)) state_next = ST_COMPUTE;
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
    // Eight baseline cores
    //-------------------------------------------------------------------------
    mm1_cores #(
        .N(BLOCK_N), .W(DATA_W), .ACC_W(ACC_W)
    ) u_mm1_cores (
        .clk(clk), .rst_n(rst_n),
        .load_en(pre_b_valid), .load_row(b_cnt[ROW_IDX_W-1:0]),
        .a11(a11), .a12(a12), .a21(a21), .a22(a22),
        .b11(b11), .b12(b12), .b21(b21), .b22(b22),
        .t_valid(t_valid),
        .q_valid(q_valid),
        .q1(q1), .q2(q2), .q3(q3), .q4(q4),
        .q5(q5), .q6(q6), .q7(q7), .q8(q8)
    );

    //-------------------------------------------------------------------------
    // Conventional post-addition (2 pipeline stages, mirrors strassen_post)
    //
    //   C11 = Q1 + Q3    C12 = Q2 + Q4    C21 = Q5 + Q7    C22 = Q6 + Q8
    //-------------------------------------------------------------------------
    logic signed [BLOCK_N*ACC_W-1:0] q1r, q2r, q3r, q4r, q5r, q6r, q7r, q8r;
    logic v1, c_valid_post;

    function automatic logic signed [CW-1:0] sextq(input logic signed [ACC_W-1:0] x);
        sextq = $signed({{(CW-ACC_W){x[ACC_W-1]}}, x});
    endfunction

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            v1 <= 1'b0; c_valid_post <= 1'b0;
            q1r <= '0; q2r <= '0; q3r <= '0; q4r <= '0;
            q5r <= '0; q6r <= '0; q7r <= '0; q8r <= '0;
        end else begin
            v1 <= q_valid;
            c_valid_post <= v1;
            q1r <= q1; q2r <= q2; q3r <= q3; q4r <= q4;
            q5r <= q5; q6r <= q6; q7r <= q7; q8r <= q8;
        end
    end

    logic signed [BLOCK_N*CW-1:0] c11_row, c12_row, c21_row, c22_row;

    always_ff @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            c11_row <= '0; c12_row <= '0; c21_row <= '0; c22_row <= '0;
        end else begin
            for (int i = 0; i < BLOCK_N; i++) begin
                c11_row[i*CW +: CW] <= sextq(q1r[i*ACC_W +: ACC_W]) + sextq(q3r[i*ACC_W +: ACC_W]);
                c12_row[i*CW +: CW] <= sextq(q2r[i*ACC_W +: ACC_W]) + sextq(q4r[i*ACC_W +: ACC_W]);
                c21_row[i*CW +: CW] <= sextq(q5r[i*ACC_W +: ACC_W]) + sextq(q7r[i*ACC_W +: ACC_W]);
                c22_row[i*CW +: CW] <= sextq(q6r[i*ACC_W +: ACC_W]) + sextq(q8r[i*ACC_W +: ACC_W]);
            end
        end
    end

    assign c_valid = c_valid_post & (state == ST_COMPUTE);
    assign c_row   = {c22_row, c21_row, c12_row, c11_row};

endmodule
