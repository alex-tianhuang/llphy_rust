This README was written by Claude Sonnet 4.6 because I hate writing documentation.
If you find a mistake, please tell me about it!

# llphy_rust

A Rust-based implementation of [LLPhyScore](https://github.com/julie-forman-kay-lab/LLPhyScore), a tool for computing phase separation propensity and related biophysical features of protein sequences. Originally developed by Hao Cai (@haocai1992) in Julie Forman-Kay's lab.

Given one or more amino acid sequences, LLPhyScore computes a set of biophysical features and returns scores as raw values, z-scores, or percentiles relative to the human IDRome.

---

## Repository Structure

```
pkg_data/              # Model data files and human IDRome scores, required by packages
rust/
  llphyscore_core/     # Core library: featurization, datatypes, scoring
  llphyscore/          # CLI binary: reads FASTA, writes CSV
  llphyscore_pylib/    # Python bindings (PyO3 + maturin)
  llphyscore_buildpd/  # Build-time tool for generating package data
python/
  llphyscore/          # Python package source
```

---

## Installation: Python

**Requirements:** Python 3.8+, pip, Rust toolchain (install from [rustup.rs](https://rustup.rs))

```bash
# Download this repo
git clone https://github.com/alex-tianhuang/llphy_rust.git
cd llphy_rust/rust/llphyscore_pylib

# Install maturin
pip install maturin

# Build and install the Python package
maturin develop --release

# Use in python
python
```

### Usage

```python
from llphyscore import LLPhyScoreCalculator

calc = LLPhyScoreCalculator()

# From a list of sequences
results = calc.calculate(["MSEQVENIIDLKKELQKTLTEEFKQVFRNLNEQQRNLIGEMLQSF"])
# Returns: list of {feature_name: score} dicts

# From a dict of named sequences
results = calc.calculate({"my_protein": "MSEQVENIIDLKKELQKTLTEEFKQVFRNLNEQQRNLIGEMLQSF"})
# Returns: {name: {feature_name: score}}

# Options:
# score_type = "raw" | "z-score" | "percentile"  (default: "percentile")
# model_training_base = "human" | "human+PDB" | "PDB"  (default: "human+PDB")
# disable_pbar = True  (suppress progress bars)
results = calc.calculate(sequences, score_type="z-score", disable_pbar=True)
```

---

## Installation: Rust (CLI)

**Requirements:** Rust toolchain (install from [rustup.rs](https://rustup.rs))

```bash
# Download this repo
git clone https://github.com/alex-tianhuang/llphy_rust.git
cd llphy_rust/rust/llphyscore

# Build using cargo, compiled executable at `../target/release/llphyscore`
cargo build --release
```

### Making `llphyscore` accessible system-wide (recommended)

If you use conda or mamba, you can copy the compiled executable directly into your active environment's `bin/` directory — no PATH changes needed, and uninstalling is as simple as removing the file:

```bash
# Activate the environment you want to install into, then:
cp ../target/release/llphyscore $CONDA_PREFIX/bin/llphyscore
```

To uninstall:
```bash
rm $CONDA_PREFIX/bin/llphyscore
```

### Usage

```
llphyscore -i <input.fasta> [-o <output.csv>] [OPTIONS]

Options:
  -i, --input-file <FILE>       Input file in FASTA format
  -o, --output-file <FILE>      Output CSV file (defaults to stdout)
  -s, --score-type <TYPE>       raw | z-score | percentile  [default: percentile]
  -m, --model-train-base <BASE> human | human+PDB | PDB     [default: human+PDB]
      --log-seq-errs-to <FILE>  Write sequence validation errors to file
      --disable-pbar            Disable progress bar
  -q, --quiet                   Suppress progress bar and sequence errors
```

---

## Models

Three models are available, each trained on a different negative dataset:

| `model_training_base` | Negative set used during training |
|---|---|
| `human` | Human proteome sequences |
| `human+PDB` | Human proteome + PDB structures (default) |
| `PDB` | PDB structures only |