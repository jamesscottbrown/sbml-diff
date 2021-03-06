# sbml-diff ([project homepage](http://sysos.eng.ox.ac.uk/tebio/sbml-diff))

sbml-diff is a tool for visually representing a single SBML model, or for visually comparing 2 or more SBML models. It can be used in three ways:

* As a command-line tool (`sbml-diff.py`)
* As a python package (see package_demo.py for an example)
* through our [web interface](http://sbml-diff.jamesscottbrown.com/upload)

It relies on elements having `id` attributes, so does not work with [SBML Level](http://sbml.org/Documents/Specifications) 1, but supports files in Level 2 Version 1 - Level 3 Version 1 Core.

## Installation

sbml-diff is written in Python 2 (in the future, it is planned to also become compatible with Python 3), and depends on the [Beautiful Soup](https://www.crummy.com/software/BeautifulSoup/) library ([installation instructions](https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-beautiful-soup)), [with the lxml parser](https://www.crummy.com/software/BeautifulSoup/bs4/doc/#installing-a-parser).

Download or ``git clone`` the code, ``cd`` into the directory, and install using ``python setup.py install``.

This will install both the package and command-line tool.

## Commandline usage

    usage: sbml-diff.py [-h] [--params] [--kinetics] [--abstract]
                        [--ignore IGNORE] [--elide ELIDE] [--colors COLORS]
                        [--labels LABELS] [--stoich] [--outfile OUTFILE]
                        [--model MODEL] [--align] [--cartoon] [--force]
                        [--hide-params] [--hide-rules] [--complete]
                        infile [infile ...]    

        Summarise one, or compare two or more, SBML models as a network or table.
        Supports five distinct kinds of output:    

        * DOT representation of reaction network (circles representing species, squares representing reactions)
        * DOT representation of an abstraction of reaction network, showing only species (--abstract)
        * DOT representation of a cartoon view of a genetic regulatory network (--cartoon)
        * a table of parameters (--params)
        * a table of kinetic laws for each reaction (--kineticstable)    

        If one or more kinds of table are requested, DOT output is not produced.    
    

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
      --align               Treat species/reactions with different ids in
                            different models as the same if they have the same set
                            of MIRIAM annnotations
      --cartoon             Draw transcription using SBOL glyphs
      --force, -f           Draw comparison even if files are identical
      --hide-params         Hide parameters modified by rules/events
      --hide-rules          Do not show rules
      --complete            If no changes, exit quietly. Otherwise return param
                            table, kinetic table, and DOT output