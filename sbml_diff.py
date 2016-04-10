from bs4 import BeautifulSoup, NavigableString
import argparse
import sys
import os

def collate_interactions(child_classifications):
    child_classifications = set(child_classifications)

    if len(child_classifications) == 1 and "constant" in child_classifications:
        return "constant"

    if "monotonic_increasing" in child_classifications and "monotonic_decreasing" in child_classifications:
        return "mixed"

    if len(child_classifications) == 2 and "monotonic_increasing" in child_classifications:
        return "monotonic_increasing"

    if len(child_classifications) == 2 and "monotonic_decreasing" in child_classifications:
        return "monotonic_decreasing"


def invert_classification(classification):

    if classification in ["mixed", "constant"]:
        return classification

    if classification == "monotonic_increasing":
        return "monotonic_decreasing"

    if classification == "monotonic_decreasing":
        return "monotonic_increasing"

def classify_basic_interaction(operator, child_classifications):
    # if any children are mixed, so is result
    if "mixed" in child_classifications:
        return "mixed"

    # monotonic function of one input:
    if operator in ["root", "exp", "ln", "log", "floor", "ceiling", "factorial"]:
        return child_classifications[0]

    if operator == "power":
        return collate_interactions(child_classifications)
        # TODO: check sign of second argument, which is in general hard

    if operator in ["plus", "times"]:
        # plus is an N-ary function
        return collate_interactions(child_classifications)

    # minus is unary or binary
    if operator == "minus":
        # unary case
        if len(child_classifications) == 1:
            invert_classification(child_classifications[0])

        # binary case
        child_classifications[1] = invert_classification(child_classifications[1])
        return collate_interactions(child_classifications)


    if operator == "divide":
        # binary operator
        child_classifications[1] = invert_classification(child_classifications[1])
        return collate_interactions(child_classifications)




def categoriseInteraction(kineticLaw, species_id):
    for math in kineticLaw.select_one("math"):
        if isinstance(math, NavigableString): continue
        return categoriseInteractionInner(math, species_id)


def categoriseInteractionInner(expression, species_id):
    # We implicitly assume constants and powers are positive.

    # Stuff we still need to handle:
    # pi, infinity, exponential2
    # delay csymbol
    # piecewise functions

    #print expression

    if expression.name == "cn":
        return "constant"

    if expression.name == "ci":
        if expression.string.strip() == species_id:
            return "monotonic_increasing"
        else:
            return "constant"


    # from a BS4 object, find the identity of the operator and array of it children
    operator = None
    args = []
    for child in expression.children:
        if not operator:
            operator = child.name
        elif child.name == "ci":
            args.append(child)

    # classify each of the children
    child_classifications = []
    for arg in args:
        child_classifications.append( categoriseInteractionInner(arg, species_id) )

    # combine these classifications based on identity of function
    return classify_basic_interaction(operator, child_classifications)


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

    # A c element may contain a species/compartment/parameter/function/reaction identifier
    # note that we do not currently handle reaction identifiers correctly
    arrows = []
    arrow_directions = []
    for reaction in model.select_one("listOfReactions").select("reaction"):
        reaction_id = reaction.attrs["id"]
        for ci in reaction.select_one("kineticLaw").select("ci"):
            species_id = ci.string.strip()
            if species_id not in species_ids: continue

            # if not a reactant, add regulatory arrow
            reactant_list, product_list, compartment = get_reaction_details(model, reaction_id)
            if species_id in reactant_list: continue
            arrows.append(species_id + "->" + reaction_id)

            arrow_direction = categoriseInteraction(reaction.select_one("kineticLaw"), species_id)
            arrow_directions.append(arrow_direction)
    return arrows, arrow_directions


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

    return reactant_list, product_list, compartment


def diff_reactions(models, colors):
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
        reactant_list, product_list, compartment = get_reaction_details(models[model_set[0]], reaction_id)

        parent_model_index = list(reaction_status[reaction_id])[0]
        parent_model = models[parent_model_index]
        reaction_name = get_reaction_name(parent_model,reaction_id)


        # one
        if len(model_set) == 1 and len(models) > 1:
            color = assign_color(models, model_set, colors)

            for reactant in reactant_list:
                print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)

            for product in product_list:
                print '%s -> %s [color="%s"];' % (reaction_id, product, color)

            if compartment not in reaction_strings.keys():
                reaction_strings[compartment] = []
            #print '%s [shape="square", color="%s"];' % (reaction_id, color)
            reaction_strings[compartment].append('%s [shape="square", color="%s", label="%s"];' % (reaction_id, color, reaction_name))

        # all
        if len(model_set) == len(models):
            if compartment not in reaction_strings.keys():
                reaction_strings[compartment] = []
            reaction_strings[compartment].append(diff_reaction_common(models, reaction_id, colors))

        # some
        if 1 < len(model_set) < len(models):
            if compartment not in reaction_strings.keys():
                reaction_strings[compartment] = []
            reaction_strings[compartment].append('%s [shape="square", color="pink", label="%s"];' % (reaction_id, reaction_name)) # TODO: compare more like all

            #print '%s -> %s [color="pink"];' % (reactant, reaction_id)

    return reaction_strings


def diff_reaction_common(models, reaction_id, colors):
    # if a reaction is shared, we need to consider whether its products, reactants and rate law are also shared

    reactant_status = {}
    product_status = {}
    rate_law = ""

    for model_num, model in enumerate(models):
        reactants, products, compartment = get_reaction_details(model, reaction_id)

        if not rate_law:
            rate_law = model.select_one("listOfReactions").find(id=reaction_id).select_one("kineticLaw")
        if rate_law != model.select_one("listOfReactions").find(id=reaction_id).select_one("kineticLaw"):
            rate_law = "different"

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

        # one
        if len(model_set) == 1 and len(models) > 1:
            color = assign_color(models, model_set, colors)
            print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)

        # all
        if len(model_set) == len(models):
            print '%s -> %s [color="grey"];' % (reactant, reaction_id)

        # some
        if 0 < len(model_set) < len(models):
            print '%s -> %s [color="pink"];' % (reactant, reaction_id)

    # product arrows
    for product in product_status:
        model_set = list(product_status[product])

        # one
        if len(model_set) == 1 and len(models) > 1:
            color = assign_color(models, model_set, colors)
            print '%s -> %s [color="%s"];' % (reaction_id, product, color)

        # all
        if len(model_set) == len(models):
            print '%s -> %s [color="grey"];' % (reaction_id, product)

        # some
        if 0 < len(model_set) < len(models):
            print '%s -> %s [color="pink"];' % (reaction_id, product)

    # rate law
    parent_model = models[model_set[0]]
    reaction_name = get_reaction_name(parent_model,reaction_id)

    if rate_law == "different":
        return '%s [shape="square", color="black", label="%s"];' % (reaction_id, reaction_name)
    else:
        return '%s [shape="square", color="grey", label="%s"];' % (reaction_id, reaction_name)


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


def diff_compartment(compartment_id, models, colors, reaction_strings):
    # add extra flag specifying status, to set color

    print "\n"
    print "subgraph cluster_%s {" % compartment_id
    print "graph[style=dotted];"
    print 'label="%s";' % compartment_id

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
        color = assign_color(models, species_status[species], colors)

        parent_model_index = list(species_status[species])[0]
        parent_model = models[parent_model_index]
        print '"%s" [color="%s",label="%s"];' % (species, color, get_species_name(parent_model, species)) # TODO: what if species name differs between models

    # for each regulatory interaction, find set of models containing it
    arrow_status = {}
    for model_num, model in enumerate(models):
        arrows, arrow_directions = get_regulatory_arrow(model, compartment_id)
        for ind, arrow in enumerate(arrows):
            if arrow not in arrow_status.keys():
                arrow_status[arrow] = set()
            arrow_status[arrow].add(model_num)

    for arrow in arrow_status:
        color = assign_color(models, arrow_status[arrow], colors)

        if arrow_directions[ind] == "monotonic_increasing":
            arrowhead = "lvee"
        elif arrow_directions[ind] == "monotonic_decreasing":
            arrowhead = "ltee"
        else:
            arrowhead = "dot"
        print '%s [color="%s", arrowhead="%s"];' % (arrow, color, arrowhead)


    print "\n"

    print "}"


def assign_color(models, model_set, colors):
    # one
    if len(model_set) == 1 and len(models) > 1:
        model_index = list(model_set)[0]
        return colors[model_index]
    # all
    elif len(model_set) == len(models):
        return "grey"

    # some - hardcoded pink


def get_reactions(model):
    reactions = []
    for r in model.select_one("listOfReactions").select("reaction"):
        reactions.append(r.attrs["id"])
    return reactions


def diff_models(models, colors, print_param_comparison=False):
    if print_param_comparison:
        compare_params(models)

    print "\n\n"
    print "digraph comparison {"

    reaction_strings = diff_reactions(models, colors)

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
        diff_compartment(compartment_id, models, colors, reaction_strings)
        # TODO: alter color if compartment_status[compartment] does not contain all models

    print "}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce graphical representation of one or more SBML models.')
    parser.add_argument('--params', '-p', help='Also print textual comparison of params', action='store_true')
    parser.add_argument('--kineticstable', help='Print textual comparison of params', action='store_true')
    parser.add_argument('--outfile', type=argparse.FileType('w'), help="Output file")
    parser.add_argument('--colors', help="Output file")
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

    # redirect STDOUT to specified file
    if args.outfile:
        sys.stdout = args.outfile

    all_models = []
    model_names = []
    for inFile in args.infile:
        html = inFile.read()
        all_models.append(BeautifulSoup(html, 'xml'))

        file_name = os.path.basename(os.path.split(inFile.name)[1])
        model_names.append(os.path.splitext(file_name)[0])


    if args.kineticstable:
        print_rate_law_table(all_models, model_names)
    else:
        diff_models(all_models, all_colors, print_param_comparison=args.params)


