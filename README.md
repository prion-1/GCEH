# Genetic Code Exploration Helper
## Visual analysis of a DNA/RNA coding sequence composition

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prion-1/GCEH/main?urlpath=%2Fdoc%2Ftree%2Fgceh.ipynb)

This tool analyzes a (coding) nucleic acid sequence with regard to its codon usage, optionally comparing it to a reference codon table of choice.
The result is then output as a text table, listing the utilized codons with their respective amino acids as well as their respecitve relative
frequencies. In addition, a graphical representation of the analysis is plotted as a stacked bar plot, highlighting the predominant codon for
an amino acid when possible (vertical alignment of the codon label).

The GC and GC3 contents are calculated as a moving average with adjustable sliding window size.

The input accepts [a, g, c, t, u] in upper/lower case. Stop codons are removed and the input sequence is truncated at the first occuring stop.
The label size on the stacked bars can be adjusted via the label size field, 0 disables the labels.

> [!TIP]
> Use: either start binder version above or launch repo in your own notebook.

Below is an example output, analyzing the codon distribution in the human PINK1 CDS (ENST00000321556.5), as well as a comparison ofthe human PARK7
codons (ENST00000338639.10) with the _H. sapiens_ consensus codon table.

![pink1_output_demo](https://github.com/user-attachments/assets/376cdaad-3e61-4eb9-8981-158fd6e707ed)

![gceh_output_demo](https://github.com/user-attachments/assets/0be508db-962e-44bf-b4ef-ed240ec0d537)

These are GC3 graphs for human TTN (ENSG00000155657), graphing two different moving averages as well as a cumulative GC3 percentage, providing
insight into compositional uniformity:

<img alt="gc3_sliding" src="https://github.com/user-attachments/assets/53d46cde-634a-4fcc-9e0c-8d50c4131999" />

<img alt="gc3_total" src="https://github.com/user-attachments/assets/d4c592fe-e898-4c9a-9c43-1fcf2c49bd22" />
