# Genomic Graph Analyzer

C++ tool for analyzing genomic variation graphs in GFA format. Implements graph cycle detection, path finding, and k-mer frequency analysis with rolling hash optimization.

## Features

- **GFA Parser**: Reads Standard GFA 1.0 files (S and L records)
- **Cycle Detection**: Identifies and removes back edges using DFS
- **DAG Conversion**: Transforms cyclic graphs to acyclic graphs
- **K-mer Analysis**: Computes k-mer frequencies using polynomial rolling hash
- **Pattern Matching**: Efficient sequence search with hash optimization

### Quick Compilation
```bash
g++ -std=c++17 -O2 -o graph_analyzer graph_analyzer.cpp
```
## Usage
```bash
./graph_analyzer <gfa_file> [kmer_length] [top_n] [source] [target]
```

## Example
```bash
./graph_analyzer example.gfa 9 10 1 100
```
## Repository Structure
```bash
genomic-graph-analyzer/
├── src/             # Source code (main.cpp)
├── example.gfa      # Example GFA file (optional)
├── .gitignore
├── LICENSE
└── README.md

```
