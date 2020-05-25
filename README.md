# Pythgroups
A Python implementation of the Pathgroups algorithm for various comparative genomic problems, based on the work of Chunfang Zheng and David Sankoff.

## Goals
Based on the original Java implementation of Pathgroups provided by the researchers (see References), we seek to improve its usability and readability, along with expanding its functionality, in our Python implementation. Our biggest planned improvements:
- Provide a config file for the user to input the tree structure using the Newick format, along with any other data needed from the user.
- Improve the clarity of the input and output of the program
- Improve the readability of the code.
- Provide implementations of the algorithm for other comparative genomic problems, such as the quartet and aliquoting problems.

## References
Zheng, C., Sankoff, D. *On the PATHGROUPS approach to rapid small phylogeny*. BMC Bioinformatics 12, S4 (2011). https://doi.org/10.1186/1471-2105-12-S1-S4 

Zheng C., *Pathgroups, a dynamic data structure for genome reconstruction problems*. Bioinformatics, Volume 26, Issue 13, Pages 1587–1594, (2010). https://doi.org/10.1093/bioinformatics/btq255

*[Software for Small Phylogeny Construction using Pathgroups](http://216.48.92.133/Softwares/smallPhylogeny/smallPhylogeny.html)*: The original Java implementation written by Zheng and Sankoff, hosted on [Sankoff's website](http://216.48.92.133)
