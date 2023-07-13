# regionperm
Association analysis of genomic regions based on permutation tests.

The use of regionperm is to answer the question: Do the regions in set A overlap with the regions in B more than expected? 
For example, in the associated publication, we used regionperm to ask the question: Are the promoter-interacting fragments that we identified overlap with regulatory histone marks more than expected?

In addition to textual output, regionperm also plot the permutation test results.



---INPUT---

We need to have the two region sets (RS), A and B, as a bed-like file:
- RS A (option -A) is a subset of a finite set of all valid regions called universe (option -U). For example, in the associated publication, RS A is a set of non-coding promoter-interacting HindIII fragments as a subset of all non-coding HindIII fragments in the genome. 

- RS B (option -B) is the feature we are interested in. For example, in the associated publication, RS B is a set of ChIP-seq peaks of regulatory histone marks from Roadmap Epigenomics dataset.



---PERMUTATION TESTS---

The core function of regionperm is to statistically evaluate the association between RS A and RS B via permutation tests.

This is done by:
1. Evaluating the original number of overlap between RS A and RS B (Original evaluation).
2. Creating a number of simulations by selecting random members of the universe, evaluating the simulated overlap each time.
3. Compute the p-value and z-score to answer whether the original evaluation is significantly different (greater or less) than would be expected based on the random simulations.



---SIMULATION OPTIONS---

The number of simulations created by regionperm can be controlled using the -n option (Default: -n 1000; 1000 random simulations).

Each time regionperm creates a simulation, regionperm selects random members of the universe that matches either:
- The original count of RS A, or
- The total genomic length of RS A
(Default: -m count; match by the original count of RS A).



---OUTPUT---

regionperm outputs 4 files:
- perm_test_result.txt ---> This is summarises the result of the permutation test (original evaluation, p-value, z-score, etc.)
- perm_test_raw_data.txt ---> This is a record of the number of overlaps in the original evaluation and all subsequent simulations.
- perm_test_figure.png ---> This is a plot of the permutation test result.
- regionperm.log ---> This is a log of the permutation test process.



---HOW TO RUN---

regionperm is run by specifying the inputs (-A, -B, -U), the simulation options (-n, -m), and the output directory (-o).

For example, in the associated publication, regionperm was run using the following code:

python regionperm.py -A PIFs.bed -B hg38_E059_H3K4me1.narrowpeak -U total_HindIII_fragments.bed -n 1000 -m length -o E059_H3K4me1_outputdir/

