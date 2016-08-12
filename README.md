# sbml-diff

sbml-diff is a tool for visually representing a single SBML model, or for visually comparing 2 or more SBML models. It can be used in three ways:

* As a command-line tool (`sbml-diff.py`)
* As a python package (see package_demo.py for an example)
* through our web interface


## Installation

Download or ``git clone`` the code, ``cd`` into the directory, and install using ``python setup.py install``.

This will install both the package and command-line tool.

## Commandline usage

    usage: sbml-diff.py [-h] [--params] [--kinetics] [--abstract]
                        [--ignore IGNORE] [--elide ELIDE] [--colors COLORS]
                        [--labels LABELS] [--stoich] [--outfile OUTFILE]
                        [--model MODEL]
                        infile [infile ...]

    Summarise one, or compare two or more, SBML models as a network or table.
    Supports four distinct kinds of output:

    * DOT representation of reaction network (circles representing species, squares representing reactions)
    * DOT representation of an abstraction of reaction network, showing only species (--abstract)
    * a table of parameters (--params)
    * a table of kinetic laws for each reaction (--kineticstable)


    positional arguments:
      infile                List of input SBML files    

    optional arguments:
      -h, --help            show this help message and exit
      --params, -p          Print textual comparison of params
      --kinetics, -k        Print textual comparison of kineticLaws
      --abstract, -a        Rather than comparing all reactions, compare abstract
                            regulatory network
      --ignore IGNORE, -i IGNORE
                            List of species to ignore (comma-separated). Works
                            with -a only
      --elide ELIDE, -e ELIDE
                            List of species to elide (comma-separated). Works with
                            -a only
      --colors COLORS, -c COLORS
                            List of colors (comma-separated)
      --labels LABELS, -l LABELS
                            Style for reaction labels (none, name, name+rate,
                            rate)
      --stoich, -s          Also label edges with stoichiometry
      --outfile OUTFILE     Output file
      --model MODEL         Make visual elements not corresponding to the n'th
                            model invisible