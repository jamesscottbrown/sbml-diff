# sbml-diff

sbml-diff is a tool for visually representing a single SBML model, or for visually comparing 2 or more SBML models. It can be used in three ways:

* As a command-line tool (`sbml-diff.py`)
* As a python package (see package_demo.py for an example)
* through our web interface


## Installation

Download or ``git clone`` the code, ``cd`` into the directory, and install using ``python setup.py install``.

This will install both the package and command-line tool.

## Commandline usage

    usage: sbml-diff.py [-h] [--params] [--kineticstable] [--outfile OUTFILE]
                        [--colors COLORS] [--reaction_labels REACTION_LABELS]
                        [--model MODEL]
                        infile [infile ...]

    Produce graphical representation of one or more SBML models.

    positional arguments:
      infile                List of input SBML files

    optional arguments:
      -h, --help            show this help message and exit
      --params, -p          Also print textual comparison of params
      --kineticstable       Print textual comparison of params
      --outfile OUTFILE     Output file
      --colors COLORS       List of colors (comma-separated)
      --reaction_labels REACTION_LABELS
                            Style for reaction labels (none, name, name+rate,
                            rate)
      --model MODEL         Make visual elements not corresponding to the n'th
                            model invisible