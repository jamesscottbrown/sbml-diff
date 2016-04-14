from effect_direction import categorise_interaction

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


def get_reactions(model):
    reactions = []
    for r in model.select_one("listOfReactions").select("reaction"):
        reactions.append(r.attrs["id"])
    return reactions


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