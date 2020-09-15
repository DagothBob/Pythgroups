# Main

### `GenomeReconstruction.py`

This is the main file to be called in the Pythgroups program. The user specifies an algorithm to use as a command line argument 
when calling this script, along with the path to a yaml file. See the manual for details on its usage.

The main difference from the original Java programs is its streamlined usage, with a yaml file (`config.yaml` by default)
for customizing the user inputs and tweeking the program's behavior.

# Pathgroups (general)

### `ChoiceStructure.py`

A `ChoiceStructure` object represents a potential "choice" step in the process of maximizing the number of cycles in 
the construction of the three breakpoint graphs. ChoiceStructures are picked based on criteria described in the 
*Priorities* section of *Pathgroups, a dynamic data structure for genome reconstruction problems*. 

### `PGMFragment.py`

A `PGMFragment` object represents a fragment as defined in *Pathgroups, a dynamic data structure for genome reconstruction problems*, which is any set of genes 
connected by red edges in a genome. A set of fragments represents the current state of the genome reconstruction. 
The object simply contains data on each end of the fragment, as well as a method to combine two fragments.

### `PGMPath.py`

A `PGMPath` object represents a path as defined in *Pathgroups, a dynamic data structure for genome reconstruction problems*, which is generalized as "any 
connected subgraph of a breakpoint graph, namely any connected part of a cycle". Simply contains the head and tail of 
the path, along with the corresponding genomes they belong to. Its only method is connecting two paths.

### `Priority.py`

A `Priority` object represents a single priority value (all of which are defined in `SmallPhylogeny.py`). Based on 
criteria described in the *Priorities* section of *Pathgroups, a dynamic data structure for genome reconstruction problems*. It seems that ChoiceStructures are 
picked based on data in this

# SmallPhylogeny

### `MedianData.py`

A `MedianData` object represents a median genome in the tree structure, containing data on neighbouring genomes, 
their PGMPaths, and gene data. Also contains relevant methods, such as getting the distance between two genomes 
(typically between the median and its neighbours), and getting the ancestor neighbours of this median.

### `MedianIteration.py`

The `MedianIteration` class is used in the optimization step of the pathgroups algorithm.

### `PGMPathForAGenome.py`

`PGMPathForAGenome` is a data structure containing all the paths for a genome.

### `SmallPhylogeny.py`

`SmallPhylogeny` is where the bulk of the work in the SmallPhylogeny algorithm is performed. The priority system as found in the 
Pathgroups algorithm is defined in this class. The main difference between this and the Java program this was adapted 
from is that the priorities are defined automatically in a loop rather than manually hardcoded for each priority, 
as there was a clear pattern found in the priority values.

### `TreeStructure.py`

`TreeStructure` represents an undirected graph with each node being a genome with either 1 or 3 neighbours. 
The structure is set by specifying each median and its three neighbouring nodes, using `set_tree_structure()`.

Along with the tree structure itself, it also stores other related information:
- All the PGM paths for each genome 
- The node representations of each gene in the first genome, where each gene has a head and tail node
- The total number of genes to be found in each genome (they must all be the same length)

# DCJRearrangements

This program checks the minimum number of DCJ operations that need to be done to transform one genome into another. 
It should be noted that the DCJRearrangements algorithm hasn't yet been formally defined, as it's based on a program 
developed by Zheng that hasn't been rigorously tested.

### `BPGDistance.py`

`BPGDistance` is used to calculate the distance between two genomes by 
counting the number of cycles in a breakpoint graph.

### `BPGPath.py`

A `BPGPath` object is a path to be used in breakpoint graphs. Simply contains data on the head and tail of the path 
along with which genomes they come from.

### `DCJOperation.py`

A `DCJOperation` object represents a double cut-and-join operation, containing related data for the operation, 
as used in DCJRearrangements.

### `DCJRearrangements.py`

`DCJRearrangements` is where the bulk of the work for the DCJ Rearrangements algorithm is done. `get_result()` 
is called for each cycle of checks until the distance between both genomes is 0, which returns the new rearrangement 
of the genome for each step.

### `GeneNode.py`

A `GeneNode` object simply represents a gene node with relevant data, such as adjacent nodes and chromosome 
information. They are used in deciding what operation to use at any given step.

# GenomeHalving

### `GroupGraph.py`

`GroupGraph` is where the bulk of the work in the GenomeHalving algorithm is performed. The priority system as found in the 
Pathgroups algorithm is defined in this class. The main difference between this and the Java program this was adapted 
from is that the priorities are defined automatically in a loop rather than manually hardcoded for each priority, 
as there was a clear pattern found in the priority values.

# GenomeAliquoting (incomplete)

### `Aliquoting.py`

An extension of GGH for n > 2. It is currently incomplete.

# Other

### `config.yaml`

The user specifies relevant data and configuration options in this YAML file. 

### `Chromosome.py`

A `Chromosome` object represents a chromosome, comprised of many `Gene` objects.

### `Gene.py`

A `Gene` object represents a gene, comprised of 2 nodes (one a head and the other a tail) and the gene's name.

### `Genome.py`

A `Genome` object represents a genome, comprised of many `Chromosome` objects

### `GenomeInString.py`

`GenomeInString` is a data structure storing a list of chromosomes represented as strings 
(where the genes will be separated by white space)

### `InputPreprocessing.py`

User input gets preprocessed here if the `use_gene_family_parser` attribute is enabled in `config.yaml`, preparing it 
for use in the pathgroups algorithm. Namely, it parses and groups the genes from a gene family file 
into their respective genomes and chromosomes, ordering them appropriately.

### `NetworkxNode.py`

A `NetworkxNode` object represents a genome node in the tree structure, its attributes being its name, ID, 
and the IDs of its neighboring genomes.
