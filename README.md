# Strassen Multisystolic Array Hardware Architectures

This repository contains the **SystemVerilog RTL implementation, functional verification, and hardware evaluation** of Strassen Multisystolic Array (SMM) hardware architectures, developed as part of an M.Tech thesis project.

The project reproduces and validates the findings of the paper:

> **"Strassen Multisystolic Array Hardware Architectures"**  
> Trevor E. Pogue and Nicola Nicolici  
> *IEEE Transactions on Very Large Scale Integration (VLSI) Systems, 2025*

---

## 📖 Project Overview

Strassen's algorithm theoretically reduces the computational complexity of matrix multiplication from \(O(n^3)\) to approximately \(O(n^{2.807})\). However, translating these savings into hardware is challenging because additional additions, data movement, and memory overhead can offset the benefits for small matrices.

This project addresses that challenge by implementing custom, highly pipelined multisystolic array architectures that translate Strassen's reduced multiplication count into potential FPGA hardware resource savings.

Both conventional matrix multiplication (MM) and Strassen multisystolic matrix multiplication (SMM) architectures are implemented and evaluated under comparable hardware conditions.

## 🏗️ Architecture & Design Features

The repository implements the proposed Strassen architectures alongside conventional baselines for fair hardware comparison.

| Feature | Description |
|---|---|
| **Baseline Multisystolic Arrays (MM)** | Conventional blocked matrix multiplication using \(8^r\) systolic cores. |
| **Strassen Multisystolic Arrays (SMM)** | Strassen-based architectures using \(7^r\) systolic cores to reduce the number of multiplications. |
| **Recursion Level \(r=1\)** | Architecture supporting 16×16 matrix tiles. |
| **Recursion Level \(r=2\)** | Architecture supporting 32×32 matrix tiles. |
| **Double Buffering (DBUF)** | Additional B buffers in the processing elements (PEs) to help hide memory latency. |
| **Overlapped Multi-Tile GEMM** | Streaming controllers that load the next B tile into a shadow buffer while the current A tile is being multiplied. |

### Core Design Concepts

- **Systolic Array Processing:** Regular, pipelined data movement through processing elements.
- **Strassen Recursion:** Reduces the number of required matrix multiplications through recursive decomposition.
- **Hardware Resource Optimization:** Evaluates the trade-off between multiplication savings and the additional logic required for Strassen's pre-addition and post-addition operations.
- **Double-Buffered Dataflow:** Supports overlapping data loading and computation for improved hardware utilization.

## ✅ Functional Verification

All implemented architectures are verified using self-checking SystemVerilog testbenches.

### Verification Methodology

- **Identity Matrix Tests:** Validate matrix data routing, indexing, and output placement.
- **Deterministic Ramp Matrices:** Exercise arithmetic operations and Strassen pre-addition and post-addition networks.
- **Recursive Matrix Packing:** Custom packing and unpacking logic maps standard two-dimensional matrices into the hardware's recursive quadrant representation.
- **Double-Buffering Verification:** Checks buffer-bank selection and data loading behavior.
- **Overlap Verification:** GEMM testbenches monitor `active_bank` signals and timing logs to confirm that B-tile loading overlaps with A-tile computation.

The testbenches are designed to detect functional mismatches and validate the correctness of the implemented datapaths and controllers.

## ⚙️ Hardware Synthesis & Evaluation

The architectures are synthesized and implemented using **Xilinx Vivado** targeting the **Xilinx Kintex UltraScale FPGA (XCKU115-FLVA1517-1-C)**.

### Evaluation Methodology

Post-implementation utilization and timing reports are generated for the individual Matrix Multiplication Unit (MXU) top modules.

The following metrics are extracted:

| Metric | Purpose |
|---|---|
| **DSP Utilization** | Measures the number of DSP resources used by the multiplication datapath. |
| **LUT Utilization** | Measures combinational logic resource consumption. |
| **Flip-Flop Utilization** | Measures sequential logic resource consumption. |
| **Maximum Operating Frequency (\(F_{\max}\))** | Estimates achievable clock frequency using post-route timing results and worst negative slack (WNS). |

Hierarchical utilization reports are used to isolate the core datapath from wrapper logic, enabling a more meaningful comparison between MM and SMM architectures.

### Key Validations

- **\(r=1\):** Validates the theoretical DSP reduction factor of approximately **1.14×** when comparing SMM1 against MM1.
- **\(r=2\):** Validates the theoretical DSP reduction factor of approximately **1.31×** when comparing SMM2 against MM2.
- **Multiplier Compute Efficiency (MCE):** Evaluates how effectively the architectures utilize their multiplier resources relative to the theoretical efficiency limits.

The theoretical multiplier reduction factor for recursion level \(r\) is:

\[
\left(\frac{8}{7}\right)^r
\]

This corresponds to approximately **1.14×** for \(r=1\) and **1.31×** for \(r=2\).

The synthesis and implementation results are used to assess how closely the RTL implementations reproduce the resource-efficiency trends reported in the reference paper.

## 🛠️ Tools & Technologies

- **HDL:** SystemVerilog
- **Simulation & Verification:** Self-checking SystemVerilog testbenches
- **Synthesis & Implementation:** Xilinx Vivado
- **Target FPGA:** Xilinx Kintex UltraScale — XCKU115-FLVA1517-1-C
- **Architecture:** Systolic Arrays, Strassen's Matrix Multiplication, Recursive Matrix Decomposition
- **Hardware Metrics:** DSP, LUT, Flip-Flop Utilization, WNS, \(F_{\max}\), and MCE

## 🎯 Project Objectives

1. Implement conventional and Strassen-based multisystolic array architectures in synthesizable SystemVerilog.
2. Reproduce the reference paper's recursive architecture and dataflow organization.
3. Verify functional correctness using self-checking testbenches.
4. Evaluate hardware resource utilization and post-route timing on a Kintex UltraScale FPGA.
5. Compare the measured hardware results against the theoretical resource savings predicted by Strassen's algorithm.

## 📌 Project Status

**M.Tech Thesis Project — Implementation, Verification, and Hardware Evaluation**

This repository focuses on RTL-level architectural reproduction and quantitative hardware evaluation of Strassen multisystolic arrays.

## 📚 Reference

Pogue, T. E., and Nicolici, N., *"Strassen Multisystolic Array Hardware Architectures,"* IEEE Transactions on Very Large Scale Integration (VLSI) Systems, 2025.
