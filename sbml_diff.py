from bs4 import BeautifulSoup


def get_params(model):
    param_ids = []
    param_values = {}

    for param in model.select_one("listOfParameters").select("parameter"):
        param_id = param.attrs["id"]
        param_ids.append(param_id)
        param_values[param_id] = param.attrs["value"]

    return set(param_ids), param_values


def compare_params(model1, model2):
    param_ids_1, param_values_1 = get_params(model1)
    param_ids_2, param_values_2 = get_params(model2)

    a_only = param_ids_1.difference(param_ids_2)
    b_only = param_ids_2.difference(param_ids_1)
    both = param_ids_1.intersection(param_ids_2)

    if a_only:
        print "\nParams in first model only:"
        print "\n".join(a_only)

    if b_only:
        print "Params in second model only:", " ".join(a_only)

    if both:
        same = []
        different = []

        for param in both:
            if param_values_1[param] == param_values_2[param]:
                same.append(param)
            else:
                different.append(param)

        if different:
            print "\nParams with different values in both models:"
            for param in different:
                print param, param_values_1[param], param_values_2[param]

        if same:
            print "\nParams with same values in both models:"
            for param in same:
                print param, param_values_1[param]


def get_species(model, compartment_id):
    # Return array of ids of species
    # model is a BeatifulSoup object, and compartmentID is a string
    ids = []
    for s in model.select_one("listOfSpecies").select("species"):
        if s.attrs["compartment"] == compartment_id:
            ids.append(s.attrs["id"])
    return ids


def print_species(species_list, status):
    if status == "a_only":
        format_sting = '"%s" [color="red"];'
    elif status == "b_only":
        format_sting = '"%s" [color="green"];'
    else:
        format_sting = '"%s" [color="grey"];'

    for s in species_list:
        print format_sting % s


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


def diff_reactions(model1, model2):
    # NB. models do not have an associated compartment!
    reactions1 = set(get_reactions(model1))
    reactions2 = set(get_reactions(model2))

    a_only, b_only, both = categorise(reactions1, reactions2)

    for r in both:
        diff_reaction_common(model1, model2, r)

    for reaction_id in a_only:
        print_reaction_unique(model1, reaction_id, "a_only")

    for r in b_only:
        print_reaction_unique(model2, reaction_id, "b_only")


def print_reaction_unique(model, reaction_id, status):
    # handles a reaction that is present in only one

    if status == "a_only":
        color = "red"
    else:
        color = "green"

    reactant_list, product_list = get_reaction_details(model, reaction_id)
    for reactant in reactant_list:
        print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)
    for product in product_list:
        print '%s -> %s [color="%s"];' % (reaction_id, product, color)
    print '%s [shape="square", color="%s"];' % (reaction_id, color)


def diff_reaction_common(model1, model2, reaction_id):
    # This handles a reaction that is present in both
    reactant_list1, product_list1 = get_reaction_details(model1, reaction_id)
    reactant_list2, product_list2 = get_reaction_details(model2, reaction_id)

    # reactant arrows
    a_only, b_only, both = categorise(set(reactant_list1), set(reactant_list2))

    for r in a_only:
        print '%s -> %s [color="red"];' % (r, reaction_id)

    for r in b_only:
        print '%s -> %s [color="green"];' % (r, reaction_id)

    for r in both:
        print '%s -> %s [color="grey"];' % (r, reaction_id)

    # product arrows
    a_only, b_only, both = categorise(set(product_list1), set(product_list2))

    for r in a_only:
        print '%s -> %s [color="red"];' % (reaction_id, r)

    for r in b_only:
        print '%s -> %s [color="green"];' % (reaction_id, r)

    for r in both:
        print '%s -> %s [color="grey"];' % (reaction_id, r)

    # rate law
    r1 = model1.select_one("listOfReactions").find(id=reaction_id).select_one("kineticLaw")
    r2 = model2.select_one("listOfReactions").find(id=reaction_id).select_one("kineticLaw")

    if r1.contents == r2.contents:
        print '%s [shape="square", color="grey"];' % reaction_id
    else:
        print '%s [shape="square", color="black"];' % reaction_id


def diff_compartment(compartment_id, model1, model2):
    # add extra flag specifying status, to set color

    print "\n"
    print "subgraph cluster_%s {" % compartment_id
    print "graph[style=dotted];"
    print 'label="%s";' % compartment_id

    # compare species
    #
    species1 = set(get_species(model1, compartment_id))
    species2 = set(get_species(model2, compartment_id))

    a_only, b_only, both = categorise(species1, species2)

    print "\n"
    print_species(a_only, 'a_only')
    print_species(b_only, 'b_only')
    print_species(both, 'both')



    print "\n"
    #
#    for r2 in model2.select_one("listOfReactions").select("reaction"):
 #       if r.contents == r2.contents:
  #          print "A Match!"

    # print species


    print "}"

def get_reactions(model):
    reactions = []
    for r in model.select_one("listOfReactions").select("reaction"):
        reactions.append(r.attrs["id"])
    return reactions


def diff_models(model1, model2):
    compare_params(model1, model2)

    print "\n\n"
    print "digraph comparison {"

    diff_reactions(model1, model2)

    for compartment in model1.select('compartment'):
        compartment_id = compartment.attrs["id"]
        compartment2 = model2.find(id=compartment_id)

        if compartment2:
            diff_compartment(compartment_id, model1, model2)
        else:
            pass
            # compartment only in A

    for compartment in model2.select('compartment'):
        compartment_id = compartment.attrs["id"]
        compartment1 = model1.find(id=compartment_id)
        if not compartment1:
            pass
            # compartment only in B

    print "}"
