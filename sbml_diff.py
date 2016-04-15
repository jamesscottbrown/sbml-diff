from bs4 import BeautifulSoup
from accessor_functions import *
from generate_dot import *
from rate_laws import *
import argparse
import sys
import os


def print_rate_law_table(models, model_names):
    print "<table>"

    print "<thead><tr><th></th><th> %s </th></tr></thead>" % ("</th><th>".join(model_names))

    # get list of all reactions in all models
    reactions = []
    for model in models:
        reactions.extend(get_reactions(model))
    reactions = list(set(reactions))
    reactions.sort()

    print "<tbody>"
    for reaction_id in reactions:
        rates = []
        for model_num, model in enumerate(models):
            r = model.select_one("listOfReactions").find(id=reaction_id)
            if r:
                math_tag = r.select_one("kineticLaw").select_one("math")
                rates.append(unicode(math_tag))
            else:
                rates.append("-")

        print "<tr> <td>%s</td> <td>%s</td></tr>" % (reaction_id, "</td><td>".join(rates))

    print "</tbody></table>"


def compare_params(models):
    # For each param, find set of models containing it, and determine whether its value is the same across all models

    param_status = {}
    param_value = {}
    for model_num, model in enumerate(models):
        param_ids, param_values = get_params(model)

        for param_id in param_ids:

            if param_id not in param_status.keys():
                param_status[param_id] = set()
                param_value[param_id] = param_values[param_id]

            if param_value[param_id] != param_values[param_id]:
                param_value[param_id] = "different"

            param_status[param_id].add(model_num)

    # one
    print "\nParameters in a single model only:"

    for model_num, model in enumerate(models):
        for param_id in param_status:
            model_list = list(param_status[param_id])
            if len(model_list) == 1 and model_list[0] == model_num and len(models) > 1:
                print "Only in model %s: %s" % (model_num, param_id)

    # all
    print "\nParameters in all models:"

    for param_id in param_status:
        model_list = list(param_status[param_id])
        value = param_value[param_id]
        if value != "different":
            value = "same"
        if len(model_list) == len(models):
            print "In all models (with %s values): %s" % (value, param_id)

    # some
    print "\nParameters in some models:"
    for param_id in param_status:
        model_list = list(param_status[param_id])
        if value != "different":
            value = "same"

        if 1 < len(model_list) < len(models):
            print "In some models (with %s values): %s (in %s)" % (value, param_id, ', '.join(model_list))


def diff_rules(models, generate_dot):
    rule_status = {}
    for model_num, model in enumerate(models):
        rules = get_rules(model)

        for rule in rules:
            if rule not in rule_status.keys():
                rule_status[rule] = set()

            rule_status[rule].add(model_num)

    rule_strings = {}

    for rule_target in rule_status:
        model_set = list(rule_status[rule])
        inputs, output, compartment, rate_law = get_rule_details(models[model_set[0]], rule_target)

        if compartment not in rule_strings.keys():
            rule_strings[compartment] = []

        rule_string = diff_rule(models, rule_target, generate_dot)
        rule_strings[compartment].append(rule_string)

    return rule_strings


def diff_rule(models, rule_id, generate_dot):
    # if a reaction is shared, we need to consider whether its products, reactants and rate law are also shared

    # 'modifiers' appear in the math expression of a rule that sets 'target'
    # a rule has only one target, whereas reaction may have multiple products

    modifier_status = {}
    target_status = set()
    rate_laws = ""

    for model_num, model in enumerate(models):
        modifiers, target, compartment, rate_law = get_rule_details(model, rule_id)

        if not target:
            return ""

        if rate_law and not rate_laws:
            rate_laws = rate_law
        if rate_laws and rate_law and rate_laws != rate_law:
            rate_laws = "different"

        for modifier in modifiers:
            if modifier not in modifier_status.keys():
                modifier_status[modifier] = set()
            modifier_status[modifier].add(model_num)

        # for target in targets:
        target_status.add(model_num)

    # modifier arrows
    for modifier in modifier_status:
        model_set = list(modifier_status[modifier])
        generate_dot.print_rule_modifier_arrow(model_set, rule_id, modifier)

    # target arrows
    # for product in product_status:
    model_set = list(target_status)
    generate_dot.print_rule_target_arrow(model_set, rule_id, target)

    # rate law
    reaction_name = rule_id  # rules don't have name attributes

    converted_rate_law = ""
    if rate_laws and rate_laws != "different":
        converted_rate_law = convert_rate_law(rate_laws)

    return generate_dot.print_rule_node(model_set, rule_id, rate_laws, reaction_name, converted_rate_law)



def diff_reactions(models, generate_dot):
    # NB. reactions do not have an associated compartment!
    reaction_status = {}
    for model_num, model in enumerate(models):
        reactions = get_reactions(model)

        for reaction in reactions:
            if reaction not in reaction_status.keys():
                reaction_status[reaction] = set()

            reaction_status[reaction].add(model_num)

    reaction_strings = {}

    for reaction_id in reaction_status:
        model_set = list(reaction_status[reaction_id])
        reactant_list, product_list, compartment, rate_law = get_reaction_details(models[model_set[0]], reaction_id)

        if compartment not in reaction_strings.keys():
            reaction_strings[compartment] = []

        reaction_string = diff_reaction(models, reaction_id, generate_dot)
        reaction_strings[compartment].append(reaction_string)

    return reaction_strings


def diff_reaction(models, reaction_id, generate_dot):
    # if a reaction is shared, we need to consider whether its products, reactants and rate law are also shared

    reactant_status = {}
    product_status = {}
    rate_laws = ""

    for model_num, model in enumerate(models):
        reactants, products, compartment, rate_law = get_reaction_details(model, reaction_id)

        if rate_law and not rate_laws:
            rate_laws = rate_law
        if rate_laws and rate_law and rate_laws != rate_law:
            rate_laws = "different"

        for reactant in reactants:
            if reactant not in reactant_status.keys():
                reactant_status[reactant] = set()
            reactant_status[reactant].add(model_num)

        for product in products:
            if product not in product_status.keys():
                product_status[product] = set()
            product_status[product].add(model_num)

    # reactant arrows
    for reactant in reactant_status:
        model_set = list(reactant_status[reactant])
        generate_dot.print_reactant_arrow(model_set, reaction_id, reactant)

    # product arrows
    for product in product_status:
        model_set = list(product_status[product])
        generate_dot.print_product_arrow(model_set, reaction_id, product)

    # rate law
    parent_model = models[model_set[0]]
    reaction_name = get_reaction_name(parent_model, reaction_id)

    converted_rate_law = ""
    if rate_laws and rate_laws != "different":
        converted_rate_law = convert_rate_law(rate_laws)

    return generate_dot.print_reaction_node(model_set, reaction_id, rate_laws, reaction_name, converted_rate_law)


def diff_compartment(compartment_id, models, reaction_strings, rule_strings, generate_dot):
    # add extra flag specifying status, to set color
    generate_dot.print_compartment_header(compartment_id)

    # print the reaction squares that belong in this compartment
    if compartment_id in reaction_strings.keys():
        print "\n".join(reaction_strings[compartment_id])

    # print the rule nodes that belong in this compartment
    if compartment_id in rule_strings.keys():
        print "\n".join(rule_strings[compartment_id])

    # For each species, find set of models containing it
    species_status = {}
    for model_num, model in enumerate(models):
        for species in get_species(model, compartment_id):

            if species not in species_status.keys():
                species_status[species] = set()
            species_status[species].add(model_num)

    for species in species_status:
        parent_model_index = list(species_status[species])[0]
        parent_model = models[parent_model_index]
        species_name = get_species_name(parent_model, species)
        generate_dot.print_species_node(species_status[species], species, species_name)

    # for each regulatory interaction - (reactant, reaction, effect direction) tuple - find set of models containing it
    arrow_status = {}
    for model_num, model in enumerate(models):
        arrows = get_regulatory_arrow(model, compartment_id)

        for ind, arrow in enumerate(arrows):
            if arrow not in arrow_status.keys():
                arrow_status[arrow] = set()
            arrow_status[arrow].add(model_num)

    for ind, arrow in enumerate(arrow_status):
        arrow_parts = arrow.split('-')
        arrow_main = '-'.join(arrow_parts[:-1])
        arrow_direction = arrow_parts[-1]
        generate_dot.print_regulatory_arrow(arrow_status[arrow], arrow_main, arrow_direction)

    generate_dot.print_compartment_footer()


def diff_models(models, generate_dot, print_param_comparison=False):
    if print_param_comparison:
        compare_params(models)

    generate_dot.print_header()

    reaction_strings = diff_reactions(models, generate_dot)
    rule_strings = diff_rules(models, generate_dot)

    if "NONE" in reaction_strings.keys():
        print reaction_strings["NONE"]
    if "NONE" in rule_strings.keys():
        print rule_strings["NONE"]

    # For every compartment in any model, record which models contain it
    compartment_status = {}
    for model_num, model in enumerate(models):
        for compartment in model.select('compartment'):
            compartment_id = compartment.attrs["id"]

            if compartment_id not in compartment_status.keys():
                compartment_status[compartment_id] = set()

            compartment_status[compartment_id].add(model_num)

    for compartment_id in compartment_status:
        diff_compartment(compartment_id, models, reaction_strings, rule_strings, generate_dot)

    generate_dot.print_footer()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce graphical representation of one or more SBML models.')
    parser.add_argument('--params', '-p', help='Also print textual comparison of params', action='store_true')
    parser.add_argument('--kineticstable', help='Print textual comparison of params', action='store_true')
    parser.add_argument('--outfile', type=argparse.FileType('w'), help="Output file")
    parser.add_argument('--colors', help="List of colors (comma-separated)")
    parser.add_argument('--reaction_labels', help="Style for reaction labels (none, name, name+rate, rate)")
    parser.add_argument('--model', help="Make visula elements not corresponding to the n'th model invisible")
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
    if args.reaction_labels:
        reaction_labels = args.reaction_labels

    selected_model = ""
    if args.model:
        selected_model = args.model

    # redirect STDOUT to specified file
    if args.outfile:
        sys.stdout = args.outfile

    all_models = []
    all_model_names = []
    for inFile in args.infile:
        html = inFile.read()
        all_models.append(BeautifulSoup(html, 'xml'))

        file_name = os.path.basename(os.path.split(inFile.name)[1])
        all_model_names.append(os.path.splitext(file_name)[0])

    if args.kineticstable:
        print_rate_law_table(all_models, all_model_names)
    else:
        diff_models(all_models, GenerateDot(all_colors, num_files, reaction_label=reaction_labels, selected_model=selected_model), print_param_comparison=args.params)
