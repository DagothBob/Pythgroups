# Pythgroups
A Python implementation of the Pathgroups algorithm for various genome reconstruction problems, based on the work of Chunfang Zheng and David Sankoff. See the [manual](./manual.md) for how to install and use Pythgroups.

## Overview
Zheng and Sankoff developed Java implementations for 2 genome reconstruction problems using their Pathgroups algorithm: small phylogeny and guided genome halving (see References). They also wrote one more program without a formally defined algorithm which calculates what double cut-and-join (DCJ) operations need to be made between two genomes to minimize the break point graph (BPG) distance between them. 

What Pythgroups does is integrate each of these programs into one Python package while improving the usability and readability of the program. Some notable features include:

- A config file for the required inputs along with customization options
- A newick tree parser to be used for small phylogeny
- Cleaner output formats for each algorithm
- Extended customization options include:
  - displaying a matplotlib graph for the calculated distances between each genome for small phylogeny
  - displaying the calculated number of each DCJ operation performed between each genome in small phylogeny after distances are calculated
- Many other parameters can be modified to tweak the program behavior and output

## References
Zheng, C., Sankoff, D. *On the PATHGROUPS approach to rapid small phylogeny*. BMC Bioinformatics 12, S4 (2011). https://doi.org/10.1186/1471-2105-12-S1-S4 

Zheng C., *Pathgroups, a dynamic data structure for genome reconstruction problems*. Bioinformatics, Volume 26, Issue 13, Pages 1587–1594, (2010). https://doi.org/10.1093/bioinformatics/btq255

#### The original programs written by Zheng and Sankoff
- *[Software for Small Phylogeny Construction using Pathgroups](http://216.48.92.133/Softwares/smallPhylogeny/smallPhylogeny.html)*
- *[Software for Guided genome halving](http://216.48.92.133/Softwares/GuidedGenomeHalving/GGH.html)*

Each of these programs are hosted on [Sankoff's website](http://216.48.92.133)

## Copyright
Biopython is used under the Biopython Licence Agreement (c)2020
