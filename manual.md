# Installation 

Note that: 1. Python 3.6 or later is required for this program, and 2. if you're a Linux user, the following 
installation steps may vary depending on your distribution, but the general principles should still apply. 
It's also worth mentioning that Pythgroups is compatible with PyPy, which speeds up the program considerably but 
requires additional installation steps. See <https://www.pypy.org/download.html> for download and installation 
instructions (be sure to use PyPy3.6)

In this guide we recommend changing your working directory to `Pythgroups`:

	$ cd [path to Pythgroups directory]

All commands in this guide will assume this is your working directory.

### Without a virtual environment

The simplest way to use Pythgroups is to install the dependencies directly on your machine. This may mess with 
other Python projects or lead to clutter, but if that's not a problem for you, then just run the following command:

    $ pip install -r requirements.txt

All the dependencies will be installed directly on your machine, allowing you to run the program outside of a virtual 
environment. Details of program usage are described in the *Usage* section of this manual.

### With a virtual environment

If you don’t want to install all the dependencies of this project directly on your base system, you can set up a 
virtual environment to install the dependencies in and run the program from. We recommend doing this instead to avoid 
conflicts with other Python projects on your machine. 

Here's a more in depth guide on using pip and virtual environments in case our guide doesn't work for you: 
<https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>

To get started, run the following command if you don't already have `virtualenv` installed on your machine:

    $ pip install virtualenv

Now we'll create the virtual environment directory inside the Pythgroups directory using the following command:

    $ python3 -m venv env     (Mac/Linux)
    > python -m venv env        (Windows)

We recommend creating this in the Pythgroups directory, but you could technically create it anywhere you want. `env` 
is the name of the virtual environment directory you wish to create, which can be named whatever you want, but we’ll use 
`env` for the rest of this section.

Now we are ready to activate the virtual environment with the following command:

    $ source env/bin/activate     (Mac/Linux)
    > env\Scripts\activate.bat      (Windows) 

For Linux users: if there isn't a `bin` directory, look for the directory containing `activate` executables.

Note that now the virtual environment's name will show on the command line in brackets indicating that you are now 
running in the virtual environment. Now we can install Pythgroup's dependencies in the virtual environment with the 
following command:

    (env) $ pip install -r requirements.txt

To exit the virtual environment, simply use the `deactivate` command:
    
    (env) $ deactivate

Now our virtual environment is all set up! Now, whenever you want to use Pythgroups through the virtual environment, 
simply run the following command as described above before using Pythgroups:

    $ source env/bin/activate     (Mac/Linux)
    > env\Scripts\activate.bat        (Windows)


# Usage

NOTE: Pythgroups is currently in development. We will update this manual as more functionality gets added.

### yaml config file

Before running the program, you have to make sure to edit the desired yaml file (default: `config.yaml`) with the 
appropriate information. Check out `TestData/InputData` for examples of each algorithm's input, and 
`TestData/InputData/ConfigReference.txt` for the required config settings for each example.

##### General
`algorithm` 
> Which algorithm you wish to use. Algorithms include `SmallPhylogeny`, `DCJRearrangements`, and `GenomeHalving`

`genome_file` 
> The path to your genome data file, relative to the Pythgroups directory. For example: 
>`TestData/InputData/SmallPhylogeny2.txt`. See examples in the `TestData` directory for the required 
>format for each algorithm.

`use_gene_family_parser`
> An experimental feature that parses raw gene family data directly, rather than requiring the genome data to be 
> formatted for pathgroups. The parser expects the following data for each gene (in the given order):
> `geneName, geneFamilyID, chr, start, end, strand, dummy1, dummy2, CoGeID, dummy3`, with each line representing a gene.
> An example showing the format for each gene's data (note that they must be delimited by tab characters):

    PAC:17829606	2	1	10068084	10075187	-	25	-1	19990	1

##### SmallPhylogeny
`tree_structure` 
> The tree structure in Newick format to be used in the SmallPhylogeny algorithm. For example: 
>`(ancestor(B.rapa,B.oleracea,B.nigra))`.
> Note that the non-median genome names must match those found in the genome file, while the median genomes can be any 
> name that isn't blank

`show_diagram` 
> Whether to print out a matplotlib diagram of the tree structure along with the distances between each 
>node at the end of the program

`show_DCJR`
> Whether to run the DCJRearrangements algorithm for each node adjacency alongside their distances

`optimization_rounds`
> The amount of optimization rounds to perform. The higher the number, the more accurate but takes longer to compute.

##### DCJRearrangements
`operations`  
> List of operations available to use in the DCJRearrangements algorithm. Note that this algorithm sometimes gives 
> errors when fission and fusion operations are enabled. The following example uses only inversions and translocations: 
```
    - inversion
    - translocation
    # - fission
    # - fusion
```
`verbose_output`
> Whether to print out its calculations at each step of the DCJRearrangements algorithm, or just print out a 
>summary of the result at the end.

`minimum_chromosome`  
> The minimum number of chromosomes to maintain for fusion operations. Must be at least 1. Default: 1

`maximum_chromosome`  
> The maximum number of chromosomes for fission operations to produce. Must be at least 2. Default: 30

`which_chromosome`
> Which chromosome to operate on: -1 for random, -2 for all of them. Default: -2

`number_of_operations`  
> The maximum number of operations to perform. Must be at least 1. Default: 10.

##### GenomeHalving
`genome_to_replace`  
> Which genome to replace. `1` to replace the first genome, `2` to replace the second genome, and `0` to pick a 
>random one to replace. Default: 1.

### GenomeReconstruction.py

`GenomeReconstruction.py` is the script used to run the Pythgroups program. Note that your working directory must be 
in `Pythgroups` for it to work. Here's how to run it:

    $ python GenomeReconstruction.py [config path]

where \[config path\] is the path of your yaml file. If left empty it will default to 'config.yaml'.

Check out `TestData/OutputData` for examples of each algorithm's output
