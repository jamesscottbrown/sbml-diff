from bs4 import BeautifulSoup
import argparse
import sys


def get_params(model):
    param_ids = []
    param_values = {}

    for param in model.select_one("listOfParameters").select("parameter"):
        param_id = param.attrs["id"]
        param_ids.append(param_id)
        param_values[param_id] = param.attrs["value"]

    return set(param_ids), param_values


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

    # print params in a single model
    print "\nParameters in a single model only:"

    for model_num, model in enumerate(models):
        for param_id in param_status:
            model_list = list(param_status[param_id])
            if len(model_list) == 1 and model_list[0] == model_num:
                print "Only in model %s: %s" % (model_num, param_id)

    # print params in all
    print "\nParameters in all models:"

    for param_id in param_status:
        model_list = list(param_status[param_id])
        value = param_value[param_id]
        if value != "different":
            value = "same"
        if len(model_list) == len(models):
            print "In all models (with %s values): %s" % (value, param_id)

    # print params in some
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


def categorise(a, b):
    a_only = a.difference(b)
    b_only = b.difference(a)
    both = a.intersection(b)

    return a_only, b_only, both


def get_reaction_details(model, reaction_id):
    reaction = model.select_one("listOfReactions").find(id=reaction_id)

    reactants = reaction.select_one("listOfReactants")
    reactant_list = []
    if reactants:
        for r in reactants.select("speciesReference"):
            reactant_list.append(r.attrs["species"])

    products = reaction.select_one("listOfProducts")
    product_list = []
    if products:
        for r in products.select("speciesReference"):
            product_list.append(r.attrs["species"])

    return reactant_list, product_list


def diff_reactions(models, colors):
    # NB. reactions do not have an associated compartment!

    reaction_status = {}
    for model_num, model in enumerate(models):
        reactions = get_reactions(model)

        for reaction in reactions:
            if reaction not in reaction_status.keys():
                reaction_status[reaction] = set()

            reaction_status[reaction].add(model_num)

    reaction_strings = []

    for reaction_id in reaction_status:
        model_set = list(reaction_status[reaction_id])

        # one
        if len(model_set) == 1:
            color = assign_color(models, model_set, colors)
            print "REACTION ID", reaction_id
            reactant_list, product_list = get_reaction_details(models[model_set[0]], reaction_id)

            for reactant in reactant_list:
                print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)

            for product in product_list:
                print '%s -> %s [color="%s"];' % (reaction_id, product, color)

            reaction_strings.append( '%s [shape="square", color="%s"];' % (reaction_id, color) )

        # all

        if len(model_set) == len(models):
            reaction_strings.append(diff_reaction_common(models, reaction_id, colors))

        # TODO: some


    # TODO: categorize by compartment
    return "\n".join(reaction_strings)


def diff_reaction_common(models, reaction_id, colors):

    model1 = models[0]
    model2 = models[1]

    # This handles a reaction that is present in both
    reactant_list1, product_list1 = get_reaction_details(model1, reaction_id)
    reactant_list2, product_list2 = get_reaction_details(model2, reaction_id)

    # reactant arrows
    a_only, b_only, both = categorise(set(reactant_list1), set(reactant_list2))

    for r in a_only:
        print '%s -> %s [color="%s"];' % (r, reaction_id, colors[0])

    for r in b_only:
        print '%s -> %s [color="%s"];' % (r, reaction_id, colors[1])

    for r in both:
        print '%s -> %s [color="grey"];' % (r, reaction_id)

    # product arrows
    a_only, b_only, both = categorise(set(product_list1), set(product_list2))

    for r in a_only:
        print '%s -> %s [color="%s"];' % (reaction_id, r, colors[0])

    for r in b_only:
        print '%s -> %s [color="%s"];' % (reaction_id, r, colors[1])

    for r in both:
        print '%s -> %s [color="grey"];' % (reaction_id, r)

    # rate law
    r1 = model1.select_one("listOfReactions").find(id=reaction_id).select_one("kineticLaw")
    r2 = model2.select_one("listOfReactions").find(id=reaction_id).select_one("kineticLaw")

    if r1.contents == r2.contents:
        return '%s [shape="square", color="grey"];' % reaction_id
    else:
        return '%s [shape="square", color="black"];' % reaction_id


def diff_compartment(compartment_id, models, colors):
    # add extra flag specifying status, to set color

    print "\n"
    print "subgraph cluster_%s {" % compartment_id
    print "graph[style=dotted];"
    print 'label="%s";' % compartment_id

    # For each species, find set of models containing it
    species_status = {}
    for model_num, model in enumerate(models):
        for species in get_species(model, compartment_id):

            if species not in species_status.keys():
                species_status[species] = set()

            species_status[species].add(model_num)

    for species in species_status:
        color = assign_color(models, species_status[species], colors)
        print '"%s" [color="%s"];' % (species, color)

    print "\n"

    print "}"

def assign_color(models, model_set, colors):
    if len(model_set) == 1:
        # entity occurs in only one model
        model_index = list(model_set)[0]
        return colors[model_index]
    elif len(model_set) == len(models):
        # entity occurs in all
        return "grey"

def get_reactions(model):
    reactions = []
    for r in model.select_one("listOfReactions").select("reaction"):
        reactions.append(r.attrs["id"])
    return reactions


def diff_models(models, colors):
    model1 = models[0]
    model2 = models[1]

    compare_params(models)

    print "\n\n"
    print "digraph comparison {"

    reaction_strings = diff_reactions(models, colors)
    print reaction_strings

    # For every compartment in any model, record which models contain it
    compartment_status = {}
    for model_num, model in enumerate(models):
        for compartment in model.select('compartment'):
            compartment_id = compartment.attrs["id"]

            if compartment_id not in compartment_status.keys():
                compartment_status[compartment_id] = set()

            compartment_status[compartment_id].add(model_num)

    for compartment_id in compartment_status:
        diff_compartment(compartment_id, models, colors)
        # TODO: alter color if compartment_status[compartment] does not contain all models

    print "}"


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Produce graphical representation of one or more SBML models.')
    parser.add_argument('--outfile', type=argparse.FileType('w'), help="Output file")
    parser.add_argument('--colors', help="Output file")
    parser.add_argument('infile', type=argparse.FileType('r'), nargs="+", help="List of input SBML files")

    args = parser.parse_args()

    if args.colors:
        colors = args.colors.split(",")
        num_files = len(args.infile)

        if len(colors) != num_files:
            print "Error: number of colors (%s) does not match number of input files (%s)\n" % (len(colors), num_files)
            parser.print_help()
            sys.exit(0)

    else:
        colors = ["red", "blue"] # A only, B only, ..., all

    # redirect STDOUT to specified file
    if args.outfile:
        sys.stdout = args.outfile

    models = []
    for file in args.infile:
        html = file.read()
        models.append(BeautifulSoup(html, 'xml'))

    diff_models(models, colors)
