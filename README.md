# Genetic Code Exploration Helper
## Visual analysis of a DNA/RNA coding sequence composition


GCEH is a lightweight analysis tool for exploring codon usage, GC/GC3 composition, and the overall structural features of a protein-coding nucleotide sequence ('complexity'). It provides descriptive statistics and visualizations, optionally comparing the input against a reference codon table.

---

## **Features**

### ** Codon Usage Analysis**

* Parses a DNA or RNA coding sequence (case-insensitive; accepts A/C/G/T/U).
* Removes stop codons and truncates the sequence at the first encountered stop.
* Generates a sorted codon table listing:

  * codon triplet
  * encoded amino acid
  * relative frequency within the sequence

### ** Visualization**

* Produces a stacked bar plot of codon usage per amino acid.
* Highlights the predominant codon for each amino acid using labels.
* Label size is adjustable; set to `0` to hide labels entirely.

### ** Complexity Analysis Module**
* The "k-mer k" field sets the k-mer size used for generating the raw complexity profile.
* "Smooth" parameter applies an adjustable sliding window to smooth the resulting curve.
* The GC-balance α slider lets you control how strongly GC content influences the complexity score (from no GC contribution at zero to GC-weighted behavior).

_normalize_sequence(seq): uppercases, converts U→T, and replaces non-ACGT with N.
_shannon_entropy_from_counts(counts): mono-nucleotide Shannon entropy (bits) over A/C/G/T.
_gc_balance(counts): GC balance score peaking at 50% GC, 0 at extremes.
_longest_run(seq): length of the longest homopolymer run (A/C/G/T only).
_kmer_richness(seq, k): for each window, count how many distinct k-length words appear (ignoring any with N) and divide by the maximum possible count, min(4^k, number of valid k-mers). The score ranges from 0 (no diversity) to 1 (all possible k-mers observed).
_periodicity_penalty(seq, max_shift): scan fixed offsets (shifts 1…max_shift), compute the fraction of positions where seq[i] equals seq[i+shift] (ignoring N/non-ACGT), and take the highest fraction as the penalty. Scores range 0–1; higher means more single-base periodicity/repeat signal.
window_metrics(seq, k, max_n_frac, gcbal_alpha): per-window composite complexity score and components (entropy, GC balance, k-mer richness, homopolymer penalty, periodicity penalty, N fraction).
compute_complexity_track(sequence, window_size, step, k, max_n_frac, gcbal_alpha, return_components, smooth): slides across a sequence to produce arrays of window start/end/mid coordinates and complexity scores (optionally smoothed) plus component tracks for plotting or downstream analysis.

### ** GC & GC3 Composition**

* Computes GC and GC3 content as a sliding-window moving average.
* Window size is user-configurable.
* Cumulative GC3 plotting to assess compositional uniformity.

### ** Reference Codon Table Comparison**

* Compare your sequence against any consensus codon-usage table (included: H. sapiens, E. coli).
* Highlights dominant codons in CDS.

---

## **Usage**

### **Option 1 — Run via Binder**

Click to launch an executable notebook environment:

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prion-1/GCEH/main?urlpath=%2Fdoc%2Ftree%2Fgceh.ipynb)

### **Option 2 — Local installation**

```bash
git clone https://github.com/prion-1/GCEH.git
cd GCEH
pip install -r requirements.txt
```

Then run the Jupyter notebook:

```bash
jupyter notebook gceh.ipynb
```

or use the `complexity_analysis.py` module directly in your own scripts.

---

## **Input Requirements**

* A valid coding sequence (DNA or RNA), containing letters **A, C, G, T, U**.
* Mixed case accepted; whitespace and newlines ignored.
* Sequence is automatically truncated at the first in-frame stop codon.

---

## **Repository Structure**

```
GCEH/
├── gceh.ipynb               # Main interactive notebook
├── complexity_analysis.py   # GC/GC3 sliding-window analysis module
├── requirements.txt         # Python dependencies
├── README.md                # This file
└── LICENSE                  # GPL-3.0 license
```

---

## **License**

This project is released under the **GPL-3.0 License**.
See `LICENSE` for details.

---

> [!TIP]
> Pasta.
