import sbml_diff
from sbml_diff import *
import os
import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Summarise one, or compare two or more, SBML models as a network or table.
    Supports four distinct kinds of output:

    * DOT representation of reaction network (circles representing species, squares representing reactions)
    * DOT representation of an abstraction of reaction network, showing only species (--abstract)
    * a table of parameters (--params)
    * a table of kinetic laws for each reaction (--kineticstable)
    """, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('--params', '-p', help='Print textual comparison of params', action='store_true')
    parser.add_argument('--kinetics', '-k', help='Print textual comparison of kineticLaws', action='store_true')

    parser.add_argument('--abstract', '-a', help='Rather than comparing all reactions, compare abstract regulatory '
                        'network', action='store_true')
    parser.add_argument('--ignore', '-i', help="List of species to ignore (comma-separated). Works with -a only")
    parser.add_argument('--elide', '-e', help="List of species to elide (comma-separated). Works with -a only")

    parser.add_argument('--colors', '-c', help="List of colors (comma-separated)")
    parser.add_argument('--labels', '-l', help="Style for reaction labels (none, name, name+rate, rate)")
    parser.add_argument('--stoich', '-s', help='Also label edges with stoichiometry', action='store_true')

    parser.add_argument('--outfile', type=argparse.FileType('w'), help="Output file")
    parser.add_argument('--model', help="Make visual elements not corresponding to the n'th model invisible")

    parser.add_argument('--align', help="Treat species/reactions with different ids in different models as the "
                        "same if they have the same set of MIRIAM annnotations", action='store_true')

    parser.add_argument('--cartoon', help="DRaw transcription using SBOL glyphs", action="store_true")

    parser.add_argument('infile', type=argparse.FileType('r'), nargs="+", help="List of input SBML files")

    args = parser.parse_args()

    num_files = len(args.infile)
    if args.colors:
        all_colors = args.colors.split(",")

        if len(all_colors) != num_files:
            print "Error: number of colors (%s) does not match number of input files (%s)\n" % (len(all_colors), num_files)
            parser.print_help()
            sys.exit(0)

    else:
        all_colors = ["#FF7F00",  "#32FF00", "#19B2FF", "#654CFF",  "#E51932", "#FFFF32"]

    reaction_labels = ""
    if args.labels:
        reaction_labels = args.labels

    selected_model = ""
    if args.model:
        selected_model = args.model

    align = False
    if args.align:
        align = True

    # redirect STDOUT to specified file
    if args.outfile:
        sys.stdout = args.outfile

    rankdir = "TB"
    cartoon = False
    if args.cartoon:
        rankdir = "LR"
        cartoon = True

    all_models = []
    all_model_names = []
    for inFile in args.infile:
        model_string = inFile.read()
        all_models.append(model_string)

        file_name = os.path.basename(os.path.split(inFile.name)[1])
        all_model_names.append(os.path.splitext(file_name)[0])

    output_formatter = sbml_diff.GenerateDot(all_colors, num_files, reaction_label=reaction_labels,
                                             selected_model=selected_model, show_stoichiometry=args.stoich,
                                             rankdir=rankdir)

    try:
        if args.kinetics:
            sbml_diff.print_rate_law_table(all_models, all_model_names)
        elif args.params:
            sbml_diff.compare_params(all_models, all_model_names)
        elif args.abstract:
            ignored = []
            if args.ignore:
                ignored = args.ignore.split(',')

            elided = []
            if args.elide:
                elided = args.elide.split(',')

            sbml_diff.diff_abstract_models(all_models, output_formatter, ignored_species=ignored, elided_species=elided, align=align)
        else:
            sbml_diff.diff_models(all_models, output_formatter, align=align, cartoon=cartoon, elided_species=[])

    except RuntimeError, e:
        sys.exit(e.args[0])
