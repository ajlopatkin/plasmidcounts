# Plasmid Count Analysis

![License](https://img.shields.io/badge/license-MIT-blue.svg)
![Python Version](https://img.shields.io/badge/python-3.8+-blue.svg)

This repository serves as the supplementary code base for the research paper titled: **Plasmid prevalence is independent of antibiotic resistance in environmental Enterobacteriaceae**, authored by: **Danya Gewurz ([dgewurz@gmail.com](mailto:dgewurz@gmail.com)), Suhyeon Kim ([skim328@ur.rochester.edu](mailto:skim328@ur.rochester.edu)), Abhishek Sharma ([ashar58@ur.rochester.edu](mailto:ashar58@ur.rochester.edu)), Joanna Harrison ([joannacharrison@comcast.net](mailto:joannacharrison@comcast.net)), Ivan Lee ([ilee24@u.rochester.edu](mailto:ilee24@u.rochester.edu)), Nicole Rondeau ([nicole.rondeau@nyulangone.org](mailto:nicole.rondeau@nyulangone.org)), JJ Miranda ([jmiranda@barnard.edu](mailto:jmiranda@barnard.edu)), Brian Mailloux ([bmailloux@barnard.edu](mailto:bmailloux@barnard.edu)), Kerry Hamilton ([kerry.hamilton@asu.edu](mailto:kerry.hamilton@asu.edu)), Allison J. Lopaltin ([allison.lopatkin@rochester.edu](mailto:allison.lopatkin@rochester.edu))**. The provided scripts and environment configurations allow users to replicate the plasmid count analysis performed in our study.

## Table of Contents

* [Project Overview](#project-overview)
* [Repository Structure](#repository-structure)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Usage](#usage)
  * [main.py](#mainpy)
  * [calibration\_bootstrap.py](#calibration_bootstrappypy)
* [Outputs](#outputs)
* [Sensitivity Parameter Generation](#sensitivity-parameter-generation)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)

## Project Overview

This project processes genomic FASTA files to detect and quantify plasmids. The provided scripts automate the process of running PlasmidFinder, aggregating the results for downstream analysis, and include utilities for automated parameter sweeping and performance evaluation to optimize plasmid detection parameters.

## Repository Structure

```
├── plasmidCount_env.yml
├── main.py
├── plasmidcounts.py
├── generate_sensitivity_data.py
├── calibration_bootstrap.py
├── configuration_example.yml
├── README.md
└── LICENSE
```

* **plasmidCount\_env.yml**: Conda environment file containing all dependencies required to run the scripts.
* **main.py**: The primary script to execute the plasmid count or sensitivity analysis, based on settings defined in `configuration_example.yml`.
* **plasmidcounts.py**: Core plasmid detection logic using PlasmidFinder and BLAST.
* **generate\_sensitivity\_data.py**: Executes the analysis across various parameter sets to assess parameter sensitivity.
* **calibration\_bootstrap.py**: Evaluates and visualizes performance of parameter sets using bootstrap sampling and a calibration/holdout strategy.
* **configuration\_example.yml**: Example YAML file where users specify all paths and analysis mode. This replaces manual terminal input of paths.
* **README.md**: This documentation file.
* **LICENSE**: Licensing information for the repository.

## Prerequisites

* **Operating System**: Windows, macOS, or Linux
* **Python**: Version 3.8 or higher
* **Conda**: Anaconda or Miniconda installed
* **Dataset and Database**: [Example Data and PLSDB database](https://rochester.box.com/s/u7vzsq7ov03h4zt9d4nmq2lk117qvrgc)

  * `example_data/`

    * `short_read/`: Contains 205 short-read FASTA files
    * `long_read/`: Contains 33 placeholder genomes for long-read assemblies

  * `db/` (required for BLAST-based analysis):

    * `plsdb.fna`: Main plasmid reference database
    * `plsdb.fna.nhr`, `plsdb.fna.nin`, `plsdb.fna.nsq`: BLAST-formatted index files

> Make sure `example_data/`, `db/`, `main.py`, and all scripts are located in the same directory before running the pipeline.

## Usage

### main.py

**Description**: Central interface for the entire plasmid analysis pipeline. This script no longer prompts the user for input paths via terminal. Instead, all configuration is handled through a `configuration_example.yml` file.

> **Important**: The following files **must all be located in the same folder**:
>
> * `main.py`, `plasmidcounts.py`, `generate_sensitivity_data.py`
> * `configuration_example.yml`
> * input data folders (`example_data/`, `db/`)

### Step-by-Step Instructions

0. Ensure all required files (scripts, configuration YAML, and example data folders) are in the **same folder**.

1. **Clone the Repository**

```bash
git clone https://github.com/ajlopatkin/plasmidcounts
cd plasmidcounts
```

2. **Create and Activate the Conda Environment**

```bash
conda env create -f plasmidCount_env.yml
conda activate plasmidCount_env
```

3. **Edit the Configuration File**

Open `configuration_example.yml` in any text editor and update the paths as shown below:

```yaml
mode: sensitivity                           # or 'plasmidcounts'
base_dir: /full/path/to/your/script_folder/example_data  # where input data folders are located
file_dirs: short_read,long_read             # name(s) of folders with FASTA files (comma-separated)
read_types: short,long                      # corresponding read type(s)
output_base: /full/path/to/output_root      # base output directory
output_dir: output                          # folder name to create (automatically generated if it doesn't exist)
db_path: /full/path/to/db                   # path to PLSDB and its BLAST index files
python_path: /full/path/to/python           # see below
plasmidfinder_path: /full/path/to/plasmidfinder.py  # see below
```

To find the correct paths for your environment, run the following commands in the activated conda environment:

```bash
which python
which plasmidfinder.py
```

Copy and paste the results into `python_path` and `plasmidfinder_path`, respectively.

> 📝 **Note**: The folder specified in `output_dir` will be created automatically inside `output_base` if it does not already exist.
> If `mode: sensitivity`, you can omit `output_dir`, as the script will generate 100 `iteration_*/` folders automatically.

4. **Navigate to the Script/Data Folder**

```bash
cd /path/to/folder/containing/scripts/and/data
```

5. **Run the Main Script**

```bash
python main.py
```

Once started, the script will automatically read the configuration and run either plasmidcounts or sensitivity analysis based on the selected `mode`. No further input is needed during execution.

### calibration\_bootstrap.py

**Description**: This script evaluates the outputs generated by `generate_sensitivity_data.py` to determine which parameter set most accurately matches the true number of plasmids per genome. It uses bootstrap sampling and a calibration/holdout split strategy to robustly assess performance across parameter combinations.

**Requirements**:

* Sensitivity output folders (`iteration_*/`, `base_params/`)
* A CSV file (`true_counts.csv`) with the following columns:

  * `names_unique`, `shortread_genome_group`, `shortread_genome_num`, `contains_true_value`

**Usage**:

```bash
python calibration_bootstrap.py \
  --root_dir /path/to/sensitivity_outputs \
  --true_count_file true_counts.csv \
  --extension .fna \
  --match_cutoff 0.85 \
  --num_bootstrap 100
```

## Outputs

### 1. Standard plasmidcounts (`mode: plasmidcounts`)

Saved to `output_base/output_dir/`:

* `plasmidcount_aggregated.csv`
* `plasmidcount_results.csv`
* `Plasmidfinder_results.csv`

### 2. Sensitivity analysis (`mode: sensitivity`)

Saved to `output_base/`:

* `iteration_0/` to `iteration_99/`: Each contains the 3 standard outputs plus `params.txt`
* `base_params/`: Results from the default parameter set

### 3. calibration\_bootstrap.py

Saved to current directory:

* `bootstrap_performance_plot.pdf`
* Terminal summary of top-performing condition and match %

## Sensitivity Parameter Generation

**Base Parameters:**

* `contig_len_cutoff` = 500
* `mash_cutoff` = 0.05
* `plasfinder_cov` = 80
* `plasfinder_pid` = 90
* `read_threshold` = 5
* `cover_threshold` = 1
* `blast_id_cutoff` = 90
* `blast_len_thresh` = 27000
* `overlap_cutoff_short` = 0.01
* `overlap_cutoff_long` = 0.15
* `blast_len_cutoff` = 1000
* `filter_col` = True
* `filter_std` = True
* `num_contig_thresh` = 3
* `skip_interactive` = False

**Sampling Strategy:** Parameters are randomly sampled within logical ranges (e.g., `plasfinder_pid`: 65-100, `overlap_cutoff_long`: 0.06-0.42) for 100 total parameter sets.

## Contributing

1. Fork the Repository
2. Create a Feature Branch
3. Commit and Push Changes
4. Open a Pull Request

## License

This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

For any inquiries or support, please email [allison.lopatkin@rochester.edu](mailto:allison.lopatkin@rochester.edu).

---

**Disclaimer**: This repository contains customized scripts for plasmid count analysis. Ensure that you understand the parameters and their implications before usage. For issues, please open a GitHub issue.
