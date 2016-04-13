from bs4 import BeautifulSoup, NavigableString
from effect_direction import categorise_interaction
from generate_dot import *
import argparse
import sys
import os


def get_params(model):
    param_ids = []
    param_values = {}

    for param in model.select_one("listOfParameters").select("parameter"):
        param_id = param.attrs["id"]
        param_ids.append(param_id)
        param_values[param_id] = param.attrs["value"]

    return set(param_ids), param_values


def get_regulatory_arrow(model, compartment):
    species_ids = get_species(model, compartment)

    arrows = []
    for reaction in model.select_one("listOfReactions").select("reaction"):
        reaction_id = reaction.attrs["id"]
        for ci in reaction.select_one("kineticLaw").select("ci"):

            # Check if this is a species id (it could validly be a species/compartment/parameter/function/reaction id)
            species_id = ci.string.strip()
            if species_id not in species_ids:
                continue

            # if not a reactant, add regulatory arrow
            reactant_list, product_list, compartment, rate_law = get_reaction_details(model, reaction_id)
            if species_id in reactant_list:
                continue

            arrow_direction = categorise_interaction(reaction.select_one("kineticLaw"), species_id)
            arrows.append('"%s" -> "%s" -%s' % (species_id, reaction_id, arrow_direction))

    return arrows


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


def get_species(model, compartment_id):
    # Return array of ids of species
    # model is a BeatifulSoup object, and compartmentID is a string
    ids = []
    for s in model.select_one("listOfSpecies").select("species"):
        if s.attrs["compartment"] == compartment_id:
            ids.append(s.attrs["id"])
    return ids


def get_species_compartment(model, species_id):
    species = model.select_one("listOfSpecies").find(id=species_id)
    return species.attrs["compartment"]


def get_reaction_details(model, reaction_id):
    reaction = model.select_one("listOfReactions").find(id=reaction_id)

    if not reaction:
        return [], [], False, False

    reactants = reaction.select_one("listOfReactants")
    reactant_list = []
    compartment = ""
    if reactants:
        for r in reactants.select("speciesReference"):
            species = r.attrs["species"]
            reactant_list.append(species)

            if not compartment:
                compartment = get_species_compartment(model, species)
            if compartment != get_species_compartment(model, species):
                compartment = "NONE"

    products = reaction.select_one("listOfProducts")
    product_list = []
    if products:
        for r in products.select("speciesReference"):
            species = r.attrs["species"]
            product_list.append(species)

            # if reaction has no reactants, try to categorise by products instead
            if not compartment:
                compartment = get_species_compartment(model, species)
            if compartment != get_species_compartment(model, species):
                compartment = "NONE"

    rate_law = reaction.select_one("kineticLaw")
    return reactant_list, product_list, compartment, rate_law


def convert_rate_law(kinetic_law):
    for math in kinetic_law.select_one("math"):
        if isinstance(math, NavigableString):
            continue
        return convert_rate_law_inner(math)[1]


def add_parens(term_elementary, terms):
    if not term_elementary[0]:
        terms[0] = "(%s)" % terms[0]
    if not term_elementary[1]:
        terms[1] = "(%s)" % terms[1]
    return terms


def convert_rate_law_inner(expression):

    # Stuff we still need to handle:
    # pi, infinity, exponential2
    # delay csymbol
    # piecewise functions

    elementary = False
    if expression.name in ["cn", "ci"]:
        elementary = True
        return elementary, expression.string.strip()

    # First child is operator; next are arguments
    operator = None
    args = []
    for child in expression.children:
        if not operator:
            operator = child.name
        else:
            if isinstance(child, NavigableString):
                continue
            args.append(child)

    children_converted = []
    children_elementary = []
    for arg in args:
        child_elementary, child_converted = convert_rate_law_inner(arg)
        children_converted.append(child_converted)
        children_elementary.append(child_elementary)

    if operator == "plus":
        return elementary, " + ".join(children_converted)
    elif operator == "minus":
        return elementary, " - ".join(children_converted)
    elif operator == "times":
        return elementary, " * ".join(children_converted)
    elif operator == "divide":
        children_converted = add_parens(children_elementary, children_converted)
        return elementary, "%s / %s" % (children_converted[0], children_converted[1])
    elif operator == "power":
        children_converted = add_parens(children_elementary, children_converted)
        return elementary, "%s ^ %s " % (children_converted[0], children_converted[1])
    elif operator in ["root", "exp", "ln", "log", "floor", "ceiling", "factorial"]:
        return elementary, "%s(%s)" % (operator, children_converted[0])


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

    num_models = len(models)

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
        generate_dot.print_reactant_arrow(num_models, model_set, reaction_id, reactant)

    # product arrows
    for product in product_status:
        model_set = list(product_status[product])
        generate_dot.print_product_arrow(num_models, model_set, reaction_id, product)

    # rate law
    parent_model = models[model_set[0]]
    reaction_name = get_reaction_name(parent_model, reaction_id)

    converted_rate_law = ""
    if rate_laws and rate_laws != "different":
        converted_rate_law = convert_rate_law(rate_laws)

    return generate_dot.print_reaction_node(num_models, model_set, reaction_id, rate_law, reaction_name, converted_rate_law)


def get_species_name(model, species_id):
        s = model.select_one("listOfSpecies").find(id=species_id)
        if "name" in s.attrs.keys() and s.attrs["name"]:
            return s.attrs["name"]
        else:
            return species_id


def get_reaction_name(model, reaction_id):
        r = model.select_one("listOfReactions").find(id=reaction_id)
        if "name" in r.attrs.keys() and r.attrs["name"]:
            return r.attrs["name"]
        else:
            return reaction_id


def diff_compartment(compartment_id, models, reaction_strings, generate_dot):
    # add extra flag specifying status, to set color
    num_models = len(models)
    generate_dot.print_compartment_header(compartment_id)

    # print the reaction squares that belong in this compartment
    print "\n".join(reaction_strings[compartment_id])

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
        generate_dot.print_species_node(num_models, species_status[species], species, species_name)

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
        generate_dot.print_regulatory_arrow(num_models, arrow_status[arrow], arrow_main, arrow_direction)

    generate_dot.print_compartment_footer()


def get_reactions(model):
    reactions = []
    for r in model.select_one("listOfReactions").select("reaction"):
        reactions.append(r.attrs["id"])
    return reactions


def diff_models(models, generate_dot, print_param_comparison=False):
    if print_param_comparison:
        compare_params(models)

    generate_dot.print_header()

    reaction_strings = diff_reactions(models, generate_dot)

    if "NONE" in reaction_strings.keys():
        print reaction_strings["NONE"]

    # For every compartment in any model, record which models contain it
    compartment_status = {}
    for model_num, model in enumerate(models):
        for compartment in model.select('compartment'):
            compartment_id = compartment.attrs["id"]

            if compartment_id not in compartment_status.keys():
                compartment_status[compartment_id] = set()

            compartment_status[compartment_id].add(model_num)

    for compartment_id in compartment_status:
        diff_compartment(compartment_id, models, reaction_strings, generate_dot)

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

    if args.colors:
        all_colors = args.colors.split(",")
        num_files = len(args.infile)

        if len(all_colors) != num_files:
            print "Error: number of colors (%s) does not match number of input files (%s)\n" % (len(all_colors), num_files)
            parser.print_help()
            sys.exit(0)

    else:
        all_colors = ["red", "blue"]  # A only, B only, ...

    reaction_labels = ""
    if args.reaction_labels:
        reaction_label = args.reaction_labels

    selected_model = ""
    if args.model != "":
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
        diff_models(all_models, GenerateDot(all_colors, reaction_label=reaction_label, selected_model=selected_model), print_param_comparison=args.params)
