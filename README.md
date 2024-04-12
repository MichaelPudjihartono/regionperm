# regionperm
<img width=50% alt="regionperm_result" src="https://github.com/MichaelPudjihartono/regionperm/assets/79574643/98c61abf-115a-4178-98d2-bc30890c9ec2">

## Association analysis of genomic regions based on permutation tests

**The use of regionperm is to answer the question: Do the regions in set A overlap with the regions in B more than expected?**

For example, in the associated publication, we used regionperm to ask the question: _Are the promoter-interacting fragments **(PIFs)** that we identified overlap with regulatory histone marks more than expected?_

In addition to textual output, regionperm also plot the permutation test results.

## Basic Usage
```
regionperm: Genomic region association analysis with permutation tests. This program allows you to assess the association between a set of genomic regions and other genomic features using permutation tests.

optional arguments:
  -h, --help            show this help message and exit
  -A A                  File containing the set of regions to randomize in BED format. File A is subset of a finite set of all valid regions in the universe file. For example, a small set of genes (file A)
                        as a subset of all genes in the genome (universe file).
  -B B                  File containing the set of genomic features that will be analyzed for their association with the regions in file A in BED format.
  -U UNIVERSE, --universe UNIVERSE
                        File containing the total valid regions from which subsequent random iterations will be sampled from in BED format. File A is a subset of the universe file
  -n NUM_ITERATIONS, --num-iterations NUM_ITERATIONS
                        Number of iterations for permutation test.
  -m {count,length}, --match-by {count,length}
                        For each iteration, do you want to randomize regions by matching the count or length of the original A region set?. Default: 'count'
  -o OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Directory to write results.
```


regionperm is run by specifying the inputs (-A, -B, -U), the simulation options (-n, -m), and the output directory (-o).

For example, in the associated publication, regionperm was run using the following code:

```
python regionperm.py -A nc_PIF.bed -B hg38_E059_H3K4me1.narrowpeak -U total_nc_genomic_fragment.bed -n 1000 -m length -o E059_H3K4me1_outputdir/
```


## Input

We need to have the two region sets (RS), A and B, as a bed-like file:
- RS A (```-A```) is a subset of a finite set of all valid regions called universe (option -U). For example, in the associated publication, RS A is a set of non-coding PIFs as a subset of all non-coding HindIII fragments in the genome. 

- RS B (```-B```) is the feature we are interested in. For example, in the associated publication, RS B is a set of ChIP-seq peaks of regulatory histone marks from Roadmap Epigenomics dataset.


## Permutation Tests


The core function of regionperm is to statistically evaluate the association between RS A and RS B via permutation tests.

This is done by:
1. Evaluating the original number of overlap between RS A and RS B (Original evaluation).
2. Creating a number of simulations by selecting random members of the universe, evaluating the simulated overlap each time.
3. Compute the p-value and z-score to answer whether the original evaluation is significantly different (greater or less) than would be expected based on the random simulations.


## Simulation Options

The number of simulations created by regionperm can be controlled using the ```-n``` option (Default: ```-n 1000```; 1000 random simulations).

Each time regionperm creates a simulation, regionperm selects random members of the universe that matches either:
- The original count of RS A (```-m count```), or
- The total genomic length of RS A (```-m length```)

(Default: ```-m count```; match by the original count of RS A).



## Output
regionperm outputs a number of files and plots.

### Result files
| File                       | Description |
|----------------------------|-------------|
|```perm_test_result.txt```  | This file summarises the result of the permutation test (original evaluation, p-value, z-score, etc.). |
|```perm_test_raw_data.txt```| This file is a record of the number of overlaps in the original evaluation and all subsequent simulations. |
|```regionperm.log```        | A log file generated during the software run. It records all operational messages, errors, and other diagnostic information, which can be crucial for troubleshooting and ensuring the software is functioning correctly. |

### Result plot
| File                       | Description |
|----------------------------|-------------|
|```perm_test_figure.png```  | A plot the permutation test result|
