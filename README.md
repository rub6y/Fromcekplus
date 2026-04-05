# Skew Symmetric Forms over F₂

A Python implementation for analyzing antisymmetric bilinear forms over the finite field F₂ (GF(2)). The code finds maximal isotropic subspaces K such that β(K, K) = 0.

## Overview

This project implements the Random Step-by-Step Method (RSBSM) for finding maximal isotropic subspaces of skew-symmetric bilinear forms over the field F₂. The theoretical background and mathematical details are described in the [Latex](./Latex) folder.

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd skewSymetricForms

# Install dependencies
pip install numpy

# Optional: for development with nix
nix-shell shell.nix
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

## Project Structure

```
.
├── MyMain.py          # Main implementation
├── AGENTS.md          # Guidelines for AI agents
├── shell.nix          # Nix development environment
├── example_beta.txt   # Example input file
├── logs/              # Output logs (input mode)
├── log.txt            # Output log (default mode)
└── Latex/
    ├── main.tex       # Mathematical description
    └── Frącek_.pdf    # Reference document
```

## Classes

- `LinearSpaceF2`: Represents a linear subspace over F₂
- `LinearMapF2`: Represents a linear map between F₂ vector spaces
- `skewSymetricFormF2`: Represents the skew-symmetric bilinear form β
- `f2_rank`: Computes matrix rank over F₂
- `f2_nullspace`: Computes null space basis over F₂

## License

MIT

## References

Mathematical details and theoretical background are available in the [Latex](./Latex) folder.