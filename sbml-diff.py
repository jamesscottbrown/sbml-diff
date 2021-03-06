import sbml_diff
from sbml_diff import *
import os
import sys
import argparse
from io import StringIO, BytesIO
import codecs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="""
    Summarise one, or compare two or more, SBML models as a network or table.
    Supports five distinct kinds of output:

    * DOT representation of reaction network (circles representing species, squares representing reactions)
    * DOT representation of an abstraction of reaction network, showing only species (--abstract)
    * DOT representation of a cartoon view of a genetic regulatory network (--cartoon)
    * a table of parameters (--params)
    * a table of kinetic laws for each reaction (--kineticstable)

    If one or more kinds of table are requested, DOT output is not produced.
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

    parser.add_argument('--cartoon', help="Draw transcription using SBOL glyphs", action="store_true")

    parser.add_argument('--force', '-f', help="Draw comparison even if files are identical", action="store_true")

    parser.add_argument('--hide-params', help="Hide parameters modified by rules/events", action="store_true")
    parser.add_argument('--hide-rules', help="Do not show rules", action="store_true")

    parser.add_argument('--sympy', help="Determine arrow directions symbolically using sympy", action="store_true")

    parser.add_argument('--complete', help="If no changes, exit quietly. Otherwise return param table, kinetic table," +
                                           " and DOT output", action="store_true")

    if (sys.version_info > (3, 0)):
        file_type = argparse.FileType('r', encoding='UTF-8')
    else:
        file_type = argparse.FileType('r')
    parser.add_argument('infile', type=file_type, nargs="+", help="List of input SBML files")

    args = parser.parse_args()

    num_files = len(args.infile)
    if args.colors:
        all_colors = args.colors.split(",")

        if len(all_colors) != num_files:
            print("Error: number of colors (%s) does not match number of input files (%s)\n" %\
                  (len(all_colors), num_files))
            parser.print_help()
            sys.exit(0)

    else:
        all_colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#ffff33", "#a65628", "#f781bf", "#999999"]

    reaction_labels = ""
    if args.labels:
        reaction_labels = args.labels

    selected_model = ""
    if args.model:
        selected_model = args.model

    align = False
    if args.align:
        align = True

    show_params = True
    if args.hide_params:
        show_params = False

    hide_rules = False
    if args.hide_rules:
        hide_rules = True

    # redirect STDOUT to specified file
    if args.outfile:
        old_stdout = args.outfile
    else:
        old_stdout = sys.stdout

    if (sys.version_info > (3, 0)):
        output_buffer = StringIO()
        f = output_buffer
    else:
        output_buffer = BytesIO()
        #f = output_buffer
        f = codecs.getwriter("utf8")(output_buffer)

    sys.stdout = f

    rankdir = "TB"
    cartoon = False
    if args.cartoon:
        rankdir = "LR"
        cartoon = True

    use_sympy = False
    if args.sympy:
        use_sympy = True

    all_models = []
    all_model_names = []
    for inFile in args.infile:
        model_string = inFile.read()
        all_models.append(model_string)

        file_name = os.path.basename(os.path.split(inFile.name)[1])
        all_model_names.append(os.path.splitext(file_name)[0])

    output_formatter = sbml_diff.GenerateDot(all_colors, num_files, reaction_label=reaction_labels,
                                             selected_model=selected_model, show_stoichiometry=args.stoich,
                                             rankdir=rankdir, model_names=all_model_names)

    sd = sbml_diff.SBMLDiff(all_models, all_model_names, output_formatter, align=align, cartoon=cartoon,
                            show_params=show_params, hide_rules=hide_rules, use_sympy=use_sympy)

    if args.complete:

        sd.print_rate_law_table()
        print("")
        sd.compare_params()
        print("")
        sd.diff_models()

        if output_formatter.differences_found:
            # print results
            sys.stdout = old_stdout
            print(f.getvalue())
        else:
            # discard results
            f.close()
            sys.stdout = old_stdout

        sys.exit()

    try:

        if args.kinetics:
            sd.print_rate_law_table()
            print("")

        if args.params:
            sd.compare_params()
            print("")

        if args.abstract:
            ignored = []
            if args.ignore:
                ignored = args.ignore.split(',')

            elided = []
            if args.elide:
                elided = args.elide.split(',')

            sd.diff_abstract_models(ignored_species=ignored, elided_species=elided)

        if not (args.kinetics or args.params or args.abstract):
            sd.diff_models()

    except RuntimeError as e:
        sys.exit(e.args[0])

    # Print results
    sys.stdout = old_stdout

    explicit_comparison = args.force or args.params or args.kinetics
    if num_files == 1 or output_formatter.differences_found or explicit_comparison:
        print(f.getvalue())
    else:
        print("No structural differences found")
