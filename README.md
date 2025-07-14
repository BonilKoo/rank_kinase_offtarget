# Kinase Off-Target Ranking Tool
This tool ranks kinase off-targets based on structural similarity to a given kinase target using KiSSim fingerprints and the KLIFS database.

It allows researchers to prioritize potential off-targets by comparing kinase binding pocket similarities, which is crucial for kinase inhibitor selectivity profiling.

## Features
- KLIFS API integration to retrieve kinase structures
- KiSSim fingerprint generation for kinase binding pockets
- Pairwise similarity analysis using nan-Euclidean distance
- Off-target ranking relative to a user-specified target kinase
- Parallelization support for faster fingerprint encoding
- Error handling for fingerprint encoding failures

## Installation
### Requirements
- Python ≥ 3.8
- Required Python packages:
    - numpy=1.24.4
    - pandas=2.0.3
    - requests=2.32.3
    - scikit-learn=1.3.2
    - tqdm=4.67.1
    - opencadd=0+untagged.20.g4867fa2.dirty
    - kissim=0+untagged.19.gf685aa5.dirty

Install dependencies via:

```bash
pip install numpy pandas requests scikit-learn tqdm opencadd kissim
```

## Usage
### Command-line Arguments
| Argument   | Description                                        | Default               |
|------------|----------------------------------------------------|-----------------------|
| `-j`       | Number of parallel jobs for KiSSim encoding        | 1                     |
| `--use_save` | Path to a precomputed kinase distance matrix (CSV) | None                  |
| `--target` | Target kinase (KLIFS name) to rank off-targets against | (Required)           |
| `--output` | Output CSV file for ranked off-targets              | ranked_offtargets.csv |

### Example Run
Run from scratch:

```bash
python kinase_offtarget_ranking.py --target CDK7
```

Use precomputed kinase distance matrix:

```bash
python kinase_offtarget_ranking.py --target CDK7 --use_save kinase_distance_matrix.csv
```

## Output
The tool outputs a CSV file with the following columns:
| Column   | Description                           |
|-----------|---------------------------------------|
| `kinase`  | Kinase name (KLIFS)                   |
| `distance`| Structural distance to target kinase |

The kinases are ranked from most similar (low distance) to least similar.

## Pipeline Overview
### 1. Download Kinase Data
- Fetches all human kinase names and structures from KLIFS API

### 2. Filter Structures
- Criteria:
    - Human species
    - DFG-in conformation
    - Resolution ≤ 3 Å
    - Quality score ≥ 6

### 3. Compute KiSSim Fingerprints
- Calculates KiSSim pocket fingerprints for the selected structures.
- Handles encoding failures with retries.

### 4. Compare Structures
- Computes nan-Euclidean distance between fingerprints.
- Converts structure-level distances to kinase-level minimal distances.

### 5. Rank Off-Targets
- For the specified `--target`, ranks kinases by structural similarity (lowest distance).

## Saved Fingerprints
To avoid recomputation, you can save and reuse the kinase distance matrix with `--use_save`.

##  Applications
- Kinase selectivity profiling
- Off-target risk assessment in drug discovery
- Kinome-wide similarity analysis

## Notes
- Requires KLIFS API access. Ensure KLIFS is accessible during execution.
- Kinase names must follow KLIFS naming conventions.

## Output Example
```python-repl
kinase,distance
CDK7,0.0
CDK2,21.85548729228192
CDK9,25.32728540306132
```

## Error Handling
- Retries fingerprint generation for structures that initially fail.
- Prints detailed logs for encoding errors.

## References
- KLIFS: https://klifs.net
```
@article{kanev2021klifs,
  title={KLIFS: an overhaul after the first 5 years of supporting kinase research},
  author={Kanev, Georgi K and de Graaf, Chris and Westerman, Bart A and de Esch, Iwan JP and Kooistra, Albert J},
  journal={Nucleic acids research},
  volume={49},
  number={D1},
  pages={D562--D569},
  year={2021},
  publisher={Oxford University Press}
}
```
- KiSSim: https://github.com/volkamerlab/kissim
```
@article{sydow2022kissim,
  title={KiSSim: predicting off-targets from structural similarities in the kinome},
  author={Sydow, Dominique and A{\ss}mann, Eva and Kooistra, Albert J and Rippmann, Friedrich and Volkamer, Andrea},
  journal={Journal of Chemical Information and Modeling},
  volume={62},
  number={10},
  pages={2600--2616},
  year={2022},
  publisher={ACS Publications}
}
```
- TeachOpenCADD: https://projects.volkamerlab.org/teachopencadd/index.html
```
@article{TeachOpenCADDKinaseEdition,
    author = {Kimber, Talia B and Sydow, Dominique and Volkamer, Andrea},
    title = {{Kinase similarity assessment pipeline for off-target prediction [v1.0]}},
    journal = {Living Journal of Computational Molecular Science},
    year = {2022},
    volume = {3},
    number = {1},
    pages = {1599},
    doi = {10.1186/s13321-019-0351-x},
}
```
