//=============================================================================
// smm1_cores_dbuf.sv
//
// Seven-core double-buffered Strassen MM array for SMM1.
//
//   Q1 = T1 * S1 ... Q7 = T7 * S7
//
// Double-buffer behavior:
//   - load_en/load_row/s1..s7 load the next S matrices into shadow buffers.
//   - swap exchanges active and shadow buffers for all seven cores.
//   - t_valid/t1..t7 stream the current T matrices through active buffers.
//
// Important:
//   - All seven cores share the same swap signal.
//   - swap only when the current tile has drained.
//   - Never assert swap and load_en on the same cycle.
//=============================================================================

(* keep_hierarchy = "yes" *)
module smm1_cores_dbuf #(
    parameter int N       = 8,     // BLOCK_N
    parameter int MW      = 17,    // T/S width
    parameter int ACC_W   = 48     // Q width
) (
    input  logic                        clk,
    input  logic                        rst_n,

    //-------------------------------------------------------------------------
    // B/S shadow-buffer load interface
    //-------------------------------------------------------------------------
    input  logic                        load_en,
    input  logic [$clog2(N)-1:0]        load_row,

    input  logic signed [N*MW-1:0]      s1,
    input  logic signed [N*MW-1:0]      s2,
    input  logic signed [N*MW-1:0]      s3,
    input  logic signed [N*MW-1:0]      s4,
    input  logic signed [N*MW-1:0]      s5,
    input  logic signed [N*MW-1:0]      s6,
    input  logic signed [N*MW-1:0]      s7,

    //-------------------------------------------------------------------------
    // Buffer swap control
    //-------------------------------------------------------------------------
    input  logic                        swap,
    output logic                        active_bank,

    //-------------------------------------------------------------------------
    // A/T streaming interface
    //-------------------------------------------------------------------------
    input  logic                        t_valid,

    input  logic signed [N*MW-1:0]      t1,
    input  logic signed [N*MW-1:0]      t2,
    input  logic signed [N*MW-1:0]      t3,
    input  logic signed [N*MW-1:0]      t4,
    input  logic signed [N*MW-1:0]      t5,
    input  logic signed [N*MW-1:0]      t6,
    input  logic signed [N*MW-1:0]      t7,

    //-------------------------------------------------------------------------
    // Q output interface
    //-------------------------------------------------------------------------
    output logic                        q_valid,

    output logic signed [N*ACC_W-1:0]   q1,
    output logic signed [N*ACC_W-1:0]   q2,
    output logic signed [N*ACC_W-1:0]   q3,
    output logic signed [N*ACC_W-1:0]   q4,
    output logic signed [N*ACC_W-1:0]   q5,
    output logic signed [N*ACC_W-1:0]   q6,
    output logic signed [N*ACC_W-1:0]   q7
);

    //-------------------------------------------------------------------------
    // Internal valid signals
    //-------------------------------------------------------------------------
    logic q1_valid;
    logic q2_valid;
    logic q3_valid;
    logic q4_valid;
    logic q5_valid;
    logic q6_valid;
    logic q7_valid;

    //-------------------------------------------------------------------------
    // Active bank status (all cores should always match)
    //-------------------------------------------------------------------------
    logic [6:0] core_active_bank;

    assign active_bank = core_active_bank[0];

    //-------------------------------------------------------------------------
    // Core 1: Q1 = T1 * S1
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q1 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s1),
        .swap(swap), .active_bank(core_active_bank[0]),
        .a_valid(t_valid), .a_row(t1),
        .c_valid(q1_valid), .c_row(q1)
    );

    //-------------------------------------------------------------------------
    // Core 2: Q2 = T2 * S2
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q2 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s2),
        .swap(swap), .active_bank(core_active_bank[1]),
        .a_valid(t_valid), .a_row(t2),
        .c_valid(q2_valid), .c_row(q2)
    );

    //-------------------------------------------------------------------------
    // Core 3: Q3 = T3 * S3
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q3 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s3),
        .swap(swap), .active_bank(core_active_bank[2]),
        .a_valid(t_valid), .a_row(t3),
        .c_valid(q3_valid), .c_row(q3)
    );

    //-------------------------------------------------------------------------
    // Core 4: Q4 = T4 * S4
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q4 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s4),
        .swap(swap), .active_bank(core_active_bank[3]),
        .a_valid(t_valid), .a_row(t4),
        .c_valid(q4_valid), .c_row(q4)
    );

    //-------------------------------------------------------------------------
    // Core 5: Q5 = T5 * S5
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q5 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s5),
        .swap(swap), .active_bank(core_active_bank[4]),
        .a_valid(t_valid), .a_row(t5),
        .c_valid(q5_valid), .c_row(q5)
    );

    //-------------------------------------------------------------------------
    // Core 6: Q6 = T6 * S6
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q6 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s6),
        .swap(swap), .active_bank(core_active_bank[5]),
        .a_valid(t_valid), .a_row(t6),
        .c_valid(q6_valid), .c_row(q6)
    );

    //-------------------------------------------------------------------------
    // Core 7: Q7 = T7 * S7
    //-------------------------------------------------------------------------
    mm_ws_systolic_dbuf #(
        .X(N), .Y(N), .A_W(MW), .B_W(MW), .ACC_W(ACC_W)
    ) u_q7 (
        .clk(clk), .rst_n(rst_n),
        .load_en(load_en), .load_row(load_row), .b_row(s7),
        .swap(swap), .active_bank(core_active_bank[6]),
        .a_valid(t_valid), .a_row(t7),
        .c_valid(q7_valid), .c_row(q7)
    );

    //-------------------------------------------------------------------------
    // All seven cores have the same latency -> Q1..Q7 aligned
    //-------------------------------------------------------------------------
    assign q_valid = q1_valid & q2_valid & q3_valid & q4_valid
                   & q5_valid & q6_valid & q7_valid;

endmodule
