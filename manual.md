# Installation 

Note that if you're a Linux user, the following installation steps may vary depending on your distribution, but the general principles should still apply

### Without a virtual environment

The simplest way to use Pythgroups is to install the dependencies directly on your machine. This may mess with other Python projects or lead to clutter, but if that's not a problem for you, then just run the following command:

    $ pip install -r Pythgroups/requirements.txt

All the dependencies will be installed directly on your machine, allowing you to run the program outside of a virtual environment. Details of program usage are described in the *Usage* section of this manual.

### With a virtual environment

If you don’t want to install all the dependencies of this project directly on your base system, you can set up a virtual environment to install the dependencies in and run the program from. We recommend doing this instead to avoid conflicts with other projects on your machine. To get started, run the following command if you don't already have `virtualenv` installed on your machine:

    $ pip install virtualenv

Now we'll create the virtual environment folder using the following command:

    $ python -m virtualenv Pythgroups/venv

We recommend creating this in the Pythgroups folder, but you could technically put it anywhere you want. `venv` is the name of the virtual environment folder you wish to create, which can be named whatever you want, but we’ll use `venv` for the rest of this section.

Now we are ready to activate the virtual environment with the following command:

    $ source Pythgroups/venv/Scripts/activate     (Mac/Linux)
    > Pythgroups\venv\Scripts\activate.bat        (Windows) 

If you're running Linux, `Scripts` may be named `bin` depending on your distribution.

Note that now the virtual environment's name will show on the command line in brackets indicating that you are now running in the virtual environment. Now we can install Pythgroup's dependencies in the virtual environment with the following command:

    (venv) $ pip install -r Pythgroups/requirements.txt

To exit the virtual environment, simply use the exit command:
    
    (venv) $ exit

Now our virtual environment is all set up! Now, whenever you want to use Pythgroups through the virtual environment, simply run the following command as described above before using Pythgroups:

    $ source Pythgroups/venv/Scripts/activate     (Mac/Linux)
    > Pythgroups\venv\Scripts\activate.bat        (Windows)


# Usage

NOTE: Pythgroups is currently in development. We will update this manual as more functionality gets added.

### config.yaml

Before running the program, you have to make sure to edit `Pythgroups/config.yaml` with the appropriate information. Here’s how it looks:

```YAML
# Tree structure in Newick format, used in SmallPhylogeny
# example: "(ancestor:0(B.rapa:26,B.oleracea:28,B.nigra:70))"
tree_structure: "[tree]"

# Path to the aligned genome data file
genome_file: "[path]"
```

`tree_structure` is the tree structure in Newick format to be used in the small phylogeny problem.

`genome_file` is the text file containing all the genome data. Note that we have yet to set up the data preprocessing needed for this program’s input, so make sure that your genome file has already gone through the gene table.

### GenomeReconstruction.py

`GenomeReconstruction.py` is the script used to run the Pythgroups program. Note that your working directory must be in `Pythgroups` for it to work. Here’s a rundown of its usage:

    $ python Pythgroups/GenomeReconstruction.py [algorithm]

Where `[algorithm]` is which algorithm to use. Possible algorithms:
* `SmallPhylogeny` - [INCOMPLETE] Used for the small phylogeny problem, which includes the median and quartet sub-problems.
* `GenomeAliquoting` - [INCOMPLETE] Used for the genome aliquoting problem.
* `DCJRearrangements` - Used for performing double-cut-and-join operations on two genomes until the first matches the second, outputting how many of each operation is performed. If more than two genomes are given in the input file, it will take the first two genomes and ignore the rest. Currently this only works for Inversion and Translocation operations, and we are working on the functionality for Fission and Fusion.
