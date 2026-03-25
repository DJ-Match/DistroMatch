# DistroMatch -  Distributed k Disjoint Matchings

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Platform: Linux](https://img.shields.io/badge/Platform-Linux-blue.svg)](https://www.linux.org/)
[![MPI](https://img.shields.io/badge/MPI-supported-green.svg)](https://www.open-mpi.org/)

A high-performance framework for distributed computation of $k$ Disjoint Matchings. This framework is optimized for efficiently processing large-scale graphs, utilizing distributed algorithms to ensure scalability and performance.

**The $k$ Disjoint Matchings Problem:**  
Given an edge-weighted graph, the goal is to identify $k$ pairwise edge-disjoint matchings that maximize the total weight of the edges in one of the $k$ matchings. This problem is particularly relevant in the context of reconfigurable optical technologies used in modern data center networks.

### Distributed $k$-Matching Algorithms

| **Algorithm**                        | **Approximation Guarantee**                          | **Key Features**                                                                 |
|--------------------------------------|-----------------------------------------------------|----------------------------------------------------------------------------------|
| **kCS** | $\frac{1}{3}$ | Uses iterative proposals to heaviest neighbors with asynchronous communication. Ensures conflict-free coloring and high-quality kk-dimensional matchings.   |
| **kMM** | $\frac{1 - \varepsilon}{3 - \varepsilon}, \varepsilon \in [0,1]$ | Employs local decisions with conflict handling. Achieves the best performance among the three algorithms. Guarantees a high-quality approximation.   |
| **MRepLM**      | $1 - \left(1 - \frac{0.5}{k}\right)^k > 0.393$      | Computes $k$ matchings successively. Simple and iterative approach. Based on the [LocalMaxMatching](https://github.com/LocalMaxMatching/LocalMaxMatching) algorithm for matching.                                                    |
| **MRepS**      | $1 - \left(1 - \frac{0.5}{k}\right)^k > 0.393$      | Computes $k$ matchings successively. Simple and iterative approach. Based on the Suitor algorithm for matching.                                      |


## Paper

If you use this software, please cite:

> Kathrin Hanauer, Sophia Heck, and Stefan Schmid. 2026. **DistroMatch:
Distributed Disjoint Weighted Matchings in Demand-Aware Reconfigurable
Optical Datacenters.** In 2026 International Conference on Supercomputing
(ICS ’26), July 06–09, 2026, Belfast, United Kingdom. ACM, New York, NY,
USA. https://doi.org/10.1145/3797905.3800555 

```bibtex
@inproceedings{distromatch,
  author    = {Kathrin Hanauer, Sophia Heck, Stefan Schmid},
  title     = {DistroMatch: Distributed Disjoint Weighted Matchings in Demand-Aware Reconfigurable Optical Datacenters},
  booktitle = {2026 International Conference on Supercomputing (ICS ’26)},
  publisher = {ACM},
  year      = {2026},
  doi       = {10.1145/3797905.3800555}
}
```


If you use MRepLM also cite the following paper, as our code uses code from [LocalMaxMatching](https://github.com/LocalMaxMatching/LocalMaxMatching):

```bibtex
@inproceedings{localmaxmatching,
  author    = {Marcel Birn and Vitaly Osipov and Peter Sanders and Christian Schulz and Nodari Sitchinava},
  title     = {Efficient Parallel and External Matching},
  booktitle = {Euro-Par 2013 Parallel Processing},
  series    = {Lecture Notes in Computer Science},
  volume    = {8097},
  pages     = {659--670},
  publisher = {Springer},
  year      = {2013},
  doi       = {10.1007/978-3-642-40047-6_66}
}
```



## Installation

To build the framework, use the provided script, which utilizes cmake for compilation:

```console
./easyCompile
```
<!-- This will compile the framework and prepare it for execution. -->



## Usage

Run the framework in a terminal using the following command:

```console
mpirun -np 4 build/Release/DistroMatch -i data/graph.mtx --algo=kCS --k=8 --partitioning=random --seed=0
```


### Command-Line Options

The framework provides a variety of running parameters:

```console
Usage: build/Release/DistroMatch [--help] -i <file> [-o <file>] --algo=VARIANT --k=<int> -p VARIANT --seed=<int> [--ccm_epsilon=<double>]
  --help                                   Print help.
  -i, --input_file=<file>                  Path to graph file.
  -o, --output_file=<file>                 Path to .sol output file.
  --algo=VARIANT                           Define the algorithm to solve kEC: [MRep|CPal|kM]
  --k=<int>                                The number of matchings to compute.
  -p, --partitioning=VARIANT               Define the partitioning used: [metis|umetis|uwmetis|ukahip|random|dwfennel|none]
  --seed=<int>                             Select a seed for the random fuction. [Default: 0]
  --ccm_epsilon=<double>                   In CCM all edges of a vertex v with weight (1-ccm_epsilon)*max_weight_v are considered candidates. Range: [0,1]. Default: 0.1

```


## Reproducibility of the Experiments

The repository provides the necessary scripts, workflows, and documentation to reproduce the results presented in our study. Below, we outline the experimental setup, preparation times, and instructions for running the experiments.

## 1. Overview of Experiments

The experiments in this repository are designed to evaluate various aspects of algorithms, including:

- **Strong Scaling**: Evaluating the scalability of algorithms on parallel architectures.
- **Comparison per $k$**: Comparing performance for different values of $k$ (number of matchings).
- **Comparison on Synthetic Graphs**: Analyzing performance on synthetic R-MAT graphs.
- **Comparison on Large Graphs**: Benchmarking performance on large-scale real-world graphs.
- **Tradeoff Between Time and Weight**: Investigating the tradeoff between execution time and solution quality.

---

### 2. Download Dependencies for Experiments

The table below summarizes the dependencies to the  required for each experiment:

| **Experiment**                  | **facebook** | **rmat** | **large** | **STK** | **DJ-Match** | **KaHIP** | **METIS** |
|----------------------------------|--------------|----------|-----------|---------|--------------|------------|-----------|
| **Strong Scaling**             | ✓ (a)        | ✓ (b)    |           |         |              |            |           |
| **Comparison per $k$**           | ✓            | ✓ (b)    |           | ✓       | ✓            |            |           |
| **Comparison on Synthetic Graphs** |              | ✓        |           | ✓       | ✓            |            |           |
| **Comparison on Large Graphs**   |              |          | ✓         | ✓       |              |            |           |
| **Tradeoff Between Time and Weight** | ✓          | ✓ (b)    |           |         |              |            |           |
| **Partitioning**                 | ✓ (a)        | ✓ (b)    |           |         |              | ✓          | ✓         |

**Notes:**
- (a) Excludes the largest instance (fbB).
- (b) Only includes R-MAT graphs with $2^{20}$ vertices.

---

#### (a) Download the Data

All datasets required for the experiments can be generated or downloaded using the following script. If you only need specific datasets, you can configure the flags in the script.

```bash
./exp/get_graphs/get_graphs.sh
```

The experiments use the following datasets:

- **Facebook - Real-World Trace Data**: Downloaded from [DC Traces](https://trace-collection.net/dc-traces/).
- **RMAT - Synthetic Data**: Generated using the R-MAT model ([Graph500](https://graph500.org/)) with the [NetworKit RmatGenerator](https://networkit.github.io/dev-docs/cpp_api/classNetworKit_1_1RmatGenerator.html).
- **LARGE - Sparse Real-World Graphs**: Sourced from the [SuiteSparse Matrix Collection](https://sparse.tamu.edu/).

---

#### (b) Download Reference Software

The following software is required for running the experiments. Ensure the correct executables are available and properly built in your environment:

- **State-of-the-art Offline Algorithms**:  
  Available in the [DJ-Match repository](https://github.com/DJ-Match/DJ-Match)  
  Path to the ex ecutable: ``../DJ-Match/build/Release/DJMatch``
- **State-of-the-art Streaming Algorithms**:  
  Provided by the authors of the paper ([arXiv link](https://arxiv.org/abs/2311.02073)). This code is referred to as **STK**.  
  Path to the ex ecutable: ``../GStream-dev/build/apps/kstmatch``

Make sure these paths are correctly set up in your environment before running the experiments.

### 3. Run the experiments

Each experiment follows a standardized workflow consisting of four tasks:

1. Generate the execution statements with the correct parameter settings and organize them in a folder structure.
2. Execute the experiment by running all the statements generated in $T_1$.
3. Extract the results and store them in a single file.
4. Visually represent the results in a plot.

For tasks, we provide shell scripts for each experiment. Below is a list of the experiments and their corresponding scripts:

---

#### **Experiments and Scripts**
- **Strong Scaling**  
  If you only want to evaluate scaling without partitioning (this is included in the partitioning experiment), use the following script:
  ```bash
  ./exp/run_scaling.sh
  ```

- **Comparison per $k$**  
  Run this script to compare results for different values of $k$:
  ```bash
  ./exp/run_comp_benchmark.sh
  ```

- **Comparison on Synthetic Graphs**  
  Use this script to evaluate performance on synthetic R-MAT graphs:
  ```bash
  ./exp/run_rmat.sh
  ```

- **Comparison on Large Graphs**  
  Use this script to compare performance on large-scale graphs:
  ```bash
  ./exp/run_comp_benchmark_large.sh
  ```

- **Tradeoff Between Time and Weight**  
  Run the following script to evaluate the tradeoff between execution time and solution weight:
  ```bash
  ./exp/run_local_max_tradeoff.sh
  ```

# Contributors
Sophia Heck - University of Vienna

Kathrin Hanauer - University of Vienna

Stefan Schmid - TU Berlin, Fraunhofer SIT