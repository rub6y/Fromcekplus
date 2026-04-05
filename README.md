# Overview of Fromcek PLUS

A Python implementation for analyzing antisymmetric bilinear forms over the finite field F₂. The code tries to finde maximal subspaces K such that β(K, K) = 0.

It is using the method described in Latex/Frącek_.pdf and it is not by any means optimal. But idea is to work on that! For now we use random choices to construct K step by step. 

## Installation

```bash
# Clone the repository
git clone https://github.com/rub6y/Fromcekplus.git
cd skewSymetricForms
```

## Usage

### Basic Usage (Random Beta)

```bash
python MyMain.py
```

### Command Line Options

| Flag | Description | Default |
|------|-------------|---------|
| `-m` | Dimension of W (vector space) | 10 |
| `-n` | Dimension of V (form space) | 2 |
| `-t, --tries` | Number of iterations | 100 |
| `-d, --depth` | Search depth (random control limit) | 10 |
| `-i, --input` | Input file with beta definitions | None |

### Examples

```bash
# Basic usage with defaults
python MyMain.py

# Custom dimensions
python MyMain.py -m 6 -n 3 -t 50

# With input file (see Input File Format below)
python MyMain.py -i betas.txt -t 100 -d 20

# Show help
python MyMain.py --help
```

### Input File Format

Create a text file with beta definitions:

```
# Comments start with #
# Format: m n on first line, then n matrices (each m rows of m integers)

3 2
0 1 1
1 0 0
1 0 0
0 0 1
0 0 1
1 1 0
```

You can define multiple betas in one file - each will be processed separately.

## Output

- **Default mode**: Results saved to `log.txt`
- **Input mode**: Results saved to `logs/beta_1.txt`, `logs/beta_2.txt`, etc.

Log format:
```
========== BETA ==========
B0
<matrix>
B1
<matrix>
==========================

Max dim: 3
v1: [1 0 1 0 1 0 0 1 0 0]
v2: [0 1 1 0 0 1 1 0 1 0]
v3: [1 1 0 1 1 0 0 1 0 1]
========================================

--- #1 ---
dim: 2
v1: ...
---
```
