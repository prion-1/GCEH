# Genetic Code Exploration Helper
## Visual analysis of the codon distribution in a DNA/RNA sequence

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/prion-1/GCEH/main?urlpath=%2Fdoc%2Ftree%2Fgceh.ipynb)

This tool analyzes a (coding) nucleic acid sequence with regard to its codon usage, optionally comparing it to a reference codon table of choice.
The result is then output as a text table, listing the utilized codons with their respective amino acids as well as their respecitve relative
frequencies. In addition, a graphical representation of the analysis is plotted as a stacked bar plot, highlighting the predominant codon for
an amino acid when possible (vertical alignment of the codon label).

The input accepts [a, g, c, t, u] in upper/lower case. Stop codons are removed and the input sequence is truncated at the first occuring stop.

> [!TIP]
> Use: either start binder version above or launch __gceh.ipynb__ in your own notebook.

Below is an example output, analyzing the codon distribution in the human PINK1 CDS (ENST00000321556.5), as well as a comparison ofthe human PARK7
codons (ENST00000338639.10) with the _H. sapiens_ consensus codon table.

![pink1_output_demo](https://github.com/user-attachments/assets/d10c2a00-25fc-40a1-853d-fe0a6a4e5d99)

![gceh_output_demo](https://github.com/user-attachments/assets/0be508db-962e-44bf-b4ef-ed240ec0d537)
