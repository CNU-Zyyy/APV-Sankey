
# APV-Sankey Toolbox

## Table of Contents
- [Introduction](#introduction)
- [Requirements](#requirements)
- [Installation](#installation)
- [Functions Overview](#functions-overview)
  - [calculate_nucleotide_ratio](#calculate_nucleotide_ratio)
  - [kmer_statistics](#kmer_statistics)
  - [calculate_common_kmers](#calculate_common_kmers)
  - [calculate_common_seqs](#calculate_common_seqs)
  - [kmer_concatenation](#kmer_concatenation)
  - [cluster_tree](#cluster_tree)
  - [kmer_evolution](#kmer_evolution)
  - [kmer_mergence](#kmer_mergence)
  - [create_sankey_chart](#create_sankey_chart)
- [Output Files](#output-files)
- [Contact](#contact)
  - [Author](#Author)
  - [Email](#Email)
- [License](#license)



## Introduction
The **APV-Sankey Toolbox** is a Python-based package designed for the systematic screening and visualization of aptamer sequences. It includes a set of tools for K-mer analysis, sequence clustering, Sankey diagram generation, and more. This toolbox is tailored for researchers working with SELEX data to track aptamer evolution across multiple rounds of selection.


## Requirements
To run the APV-Sankey Toolbox, the following Python packages are required:
- `pandas`
- `numpy`
- `matplotlib`
- `seaborn`
- `scipy`
- `openpyxl`

You can install these packages using:
```bash
pip install pandas numpy matplotlib seaborn scipy openpyxl
```


## Installation
Clone this repository and navigate to the directory:
```bash
git clone <https://github.com/CNU-Zyyy/APV-Sankey>
cd APV-Sankey
```

Ensure that all required packages are installed before running any scripts.


## Functions Overview

### calculate_nucleotide_ratio
- **Description:** This function calculates the base ratio of nucleotide sequences from either `.fasta` or `.xlsx` files. It outputs the ratio of each nucleotide (A, T, C, G, and U) for each sequence.
- **Usage:**
  ```bash
  python APV-Sankey.py calculate_nucleotide_ratio --input ./data/your_data.fasta
  ```

### kmer_statistics
- **Description:** This function performs K-mer analysis on sequences, calculating the frequency of K-mers in the provided `.fasta` or `.xlsx` file.
- **Usage:**
  ```bash
  python APV-Sankey.py kmer_statistics --input ./data/your_data.fasta --output ./results/kmer_statistics.xlsx
  ```

### calculate_common_kmers
- **Description:** This function identifies common K-mers between two datasets. It outputs a list of shared K-mers along with their frequencies.
- **Usage:**
  ```bash
  python APV-Sankey.py calculate_common_kmers --common_kmers_input ./data/Kmer_data1.xlsx ./data/Kmer_data2.xlsx --common_kmers_output ./results/calculate_common_kmers.xlsx
  ```

### calculate_common_seqs
- **Description:** This function calculates the common sequences between two datasets, listing sequences that appear in both files.
- **Usage:**
  ```bash
  python APV-Sankey.py calculate_common_seqs --common_seqs_input ./data/your_data1.xlsx ./data/your_data2.xlsx --common_seqs_output ./results/calculate_common_seqs.xlsx
  ```

### kmer_concatenation
- **Description:** This function concatenates K-mers from different rounds into longer sequences, providing a holistic view of potential binding regions.
- **Usage:**
  ```bash
  python APV-Sankey.py kmer_concatenation --kmer_concatenation_input ./data/Kmer_data.xlsx --kmer_concatenation_output ./results/kmer_concatenation.xlsx
  ```

### cluster_tree
- **Description:** Constructs hierarchical clustering trees based on sequence similarities. This can be used for classifying aptamer groups.
- **Usage:**
  ```bash
  python APV-Sankey.py cluster_tree --cluster_filepath ./data/cluster_data.xlsx
  ```

### kmer_evolution
- **Description:** This function simulates the evolution of K-mers across multiple rounds of SELEX, helping to track sequence variations.
- **Usage:**
  ```bash
  python APV-Sankey.py kmer_evolution --kmer_evolution_input ./data/Kmer_data1.xlsx ./data/Kmer_data2.xlsx --kmer_evolution_output ./results/kmer_evolution.xlsx
  ```

### kmer_mergence
- **Description:** Merges overlapping K-mers from adjacent groups into unified entries and aggregates their frequencies.
- **Usage:**
  ```bash
  python APV-Sankey.py kmer_mergence --kmer_mergence_input ./data/outputs/kmer_evolution.xlsx --kmer_mergence_output ./results/kmer_mergence.xlsx
  ```

### create_sankey_chart
- **Description:** Generates Sankey diagrams to visually represent K-mer evolution and frequency distribution across rounds.
- **Usage:**
  ```bash
  python APV-Sankey.py create_sankey_chart --sankey_input ./data/Sankey.xlsx
  ```


## Output Files
- **calculate_base_ratio:** Outputs a summary of nucleotide base ratios.
- **kmer_statistics:** Outputs an Excel file listing K-mer frequencies.
- **calculate_common_kmers:** Outputs an Excel file of K-mers shared between datasets.
- **calculate_common_seqs:** Outputs an Excel file detailing common sequences and their frequencies across multiple files.
- **kmer_concatenation:** Outputs an Excel file containing the results of concatenating K-mers based on specified rules.
- **cluster_tree:** Generates a visual clustering dendrogram that represents the hierarchical clustering of sequences based on Levenshtein Edit Distances (LEDs).
- **kmer_evolution:** Outputs an Excel file summarizing the evolution of k-mers between files.
- **kmer_mergence:** Outputs an Excel file that merges k-mers from adjacent groups in the input file, including the merged k-mers and their frequencies.
- **create_sankey_chart:** Outputs a visual representation of K-mer flow across rounds.


## Contact

### Author
* ZHANG Yu, Wang Yashan, GAO Yajing, Hu Keyi, Gong Hao, JIA Haijing, ZHANG Xin*, LOU Xinhui*
* Department of Chemistry, Capital Normal University, Beijing 100048, China
* Honours College of Capital Normal University, Capital Normal University, Beijing, 100048, China

### Email
For questions or support, please contact [xinzhang@cnu.edu.cn].


## License
The source code is released under an Apache v2.0 license.