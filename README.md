# Genetic Code Exploration Helper
## Visual analysis of the codon distribution in a DNA/RNA sequence

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prion-1/GCEH/main?urlpath=%2Fdoc%2Ftree%2Fgceh.ipynb)

This tool analyzes a (coding) nucleic acid sequence with regard to its codon usage, optionally comparing it to a reference codon table of choice.
The result is then output as a text table, listing the utilized codons with their respective amino acids as well as their respecitve relative
frequencies. In addition, a graphical representation of the analysis is plotted as a stacked bar plot, highlighting the predominant codon for
an amino acid when possible.

The input accepts [a, g, c, t, u] in upper/lower case. Stop codons are removed and the input sequence is truncated at the first occuring stop.

> [!TIP]
> Use: either start binder version above or launch __gceh.ipynb__ in your own notebook.
