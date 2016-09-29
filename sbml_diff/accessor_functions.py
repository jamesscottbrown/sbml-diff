from effect_direction import categorise_interaction


def get_params(model):
    """
    Get id and value for each parameter in a single model.


    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        

    Returns
    -------
    param_ids : a set containing the id of each parameter
    param_values : a dict whose keys are the id of a parameter, and values are corresponding parameter values

    """
    param_ids = []
    param_values = {}

    if not model.select_one("listOfParameters"):
        return set(), param_values

    for param in model.select_one("listOfParameters").select("parameter"):
        param_id = param.attrs["id"]
        param_ids.append(param_id)

        param_values[param_id] = "?"
        if "value" in param.attrs.keys():
            param_values[param_id] = param.attrs["value"]

    return set(param_ids), param_values


def get_regulatory_arrow(model, compartment, reactions, species_compartments, elided_reactions=False):
    """
    Find all regulatory interactions in a particular compartment of a model, and construct an array of strings
    representing these.
    A regulatory interaction exists if a kinetic law includes a species id that is not a reactant.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        
    compartment : the id of a compartment

    elided_reactions : a list of reactions, for which the corresponding reactions are not drawn
        

    Returns
    -------
    arrows : an array, each element of which is a string representing a regulatory interaction

    """
    if not elided_reactions:
        elided_reactions = []

    species_ids = get_species(model, compartment)

    arrows = []

    if not reactions:
        return arrows

    for reaction_id in reactions:
        reaction = reactions[reaction_id]
        if reaction in elided_reactions:
            continue

        for ci in reaction.select_one("kineticLaw").select("ci"):

            # Check if this is a species id (it could validly be a species/compartment/parameter/function/reaction id)
            species_id = ci.string.strip()
            if species_id not in species_ids:
                continue

            # if not a reactant, add regulatory arrow
            reactant_list, product_list, compartment, rate_law, _, _ = get_reaction_details(model, reaction, species_compartments)
            if species_id in reactant_list:
                continue

            arrow_direction = categorise_interaction(reaction.select_one("kineticLaw"), species_id)
            arrows.append((species_id, reaction_id, arrow_direction))

    return arrows


def get_species(model, compartment_id):
    """
    Get id of all species in compartment of model.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        
    compartment_id : the id of a compartment
        

    Returns
    -------
    ids : array listing id for each species in model

    """
    ids = []
    species_list = model.select_one("listOfSpecies")
    if species_list:
        for s in species_list.select("species"):
            if s.attrs["compartment"] == compartment_id:
                ids.append(s.attrs["id"])
    return ids


def get_species_compartment(model, species_id, species_compartments):
    """
    Get the id of the compartment containing a species.
    Report params as belonging to compartment 'NONE'

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        
    species_id : id of the species
        

    Returns
    -------
    compartment : id of the compartment

    """

    if not species_compartments:
        return "NONE"

    if species_id in species_compartments.keys():
        return  species_compartments[species_id]
    else:
        return "NONE"


def get_reaction_details(model, reaction, species_compartments):
    """
    Get details of a single reaction.

    Note that reactions have no compartment, but we try to identify it based on the reactants.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model

    reaction : bs4.BeautifulSoup object of the reaction of interest


    Returns
    -------
    reactant_list : list containing id of each reaction

    product_list : list containing id of each reaction

    compartment : id of compartment

    rate_law : BautifulSoup object containing math element of kineticLaw

    reactant_stoichiometries : list of stoichiometries of reactants

    product_stoichiometries : list of stoichiometries of products

    """

    if not reaction:
        return [], [], False, False, [], []

    reactants = reaction.select_one("listOfReactants")
    reactant_list = []
    reactant_stoichiometries = []
    compartment = ""
    if reactants:
        for r in reactants.select("speciesReference"):

            if "stoichiometry" in r.attrs:
                stoich = r.attrs["stoichiometry"]
            else:
                stoich = "1"

            reactant_stoichiometries.append(stoich)

            species = r.attrs["species"]
            reactant_list.append(species)

            if not compartment:
                compartment = get_species_compartment(model, species, species_compartments)
            if compartment != get_species_compartment(model, species, species_compartments):
                compartment = "NONE"

    products = reaction.select_one("listOfProducts")
    product_list = []
    product_stoichiometries = []
    if products:
        for r in products.select("speciesReference"):

            if "stoichiometry" in r.attrs:
                stoich = r.attrs["stoichiometry"]
            else:
                stoich = "1"

            product_stoichiometries.append(stoich)

            species = r.attrs["species"]
            product_list.append(species)

            # if reaction has no reactants, try to categorise by products instead
            if not compartment:
                compartment = get_species_compartment(model, species, species_compartments)
            if compartment != get_species_compartment(model, species, species_compartments):
                compartment = "NONE"

    rate_law = reaction.select_one("kineticLaw").select_one("math")
    return reactant_list, product_list, compartment, rate_law, reactant_stoichiometries, product_stoichiometries


def get_reactions(model):
    """
    Get list containing id for every reaction in model.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        

    Returns
    -------
    reactions : list containing id for every reaction in model

    """
    reactions = []

    reaction_list = model.select_one("listOfReactions")
    if not reaction_list:
        return reactions

    for r in reaction_list.select("reaction"):
        reactions.append(r.attrs["id"])
    return reactions


def get_rule_details(model, target_id, species_compartments, draw_modifier_params=True):
    """
    Given the id of a species affected by a rule, find details of that rule.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        
    target_id : the id of the species being affected
        

    Returns
    -------
    modifiers : list containing the id of each species in the 'math' element

    target : the id of the species affected by the rule

    compartment : the id of the compartment containing the target

    rate_law : BeautifulSoup object containing the math element for the rule

    """
    rule = model.select_one("listOfRules").find(variable=target_id)

    species_ids = []
    species_list = model.select_one("listOfSpecies")
    if species_list:
        for s in species_list.select("species"):
            species_ids.append(s.attrs["id"])

    if not rule:
        return [], False, False

    target = rule.attrs["variable"]

    # get modifier details
    modifiers = []
    for ci in rule.select("ci"):

        # Check if this is a species id (it could validly be a species/compartment/parameter/function/reaction id)
        species_id = ci.string.strip()
        if species_id not in species_ids:
            continue

        modifiers.append(species_id)

    compartment = get_species_compartment(model, target, species_compartments).strip()
    rate_law = rule.select_one("math")
    return modifiers, compartment, rate_law


def get_variables_set_by_rules(model):
    """
    Find list of species affected by a rule of type assignmentRule or rateRule.


    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        

    Returns
    -------
    list containing id of each species set by a rate or assignment rule

    """
    species = []
    rule_list = model.select_one("listOfRules")

    if not rule_list:
        return []

    for r in rule_list.select("assignmentRule"):
        species.append(r.attrs["variable"])
    for r in rule_list.select("rateRule"):
        species.append(r.attrs["variable"])
    return species


def get_species_name(model, species_id):
    """
    Return the name of the species with a given id.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        
    species_id : id of the species
        

    Returns
    -------
    name of the species, if set (otherwise returns the id)

    """
    s = model.select_one("listOfSpecies").find(id=species_id)
    if "name" in s.attrs.keys() and s.attrs["name"]:
        return s.attrs["name"]
    else:
        return species_id


def get_reaction_name(model, reaction_id):
    """
    Return the name of the reaction with a given id.

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model
        
    reaction_id : id of the reaction
        

    Returns
    -------
    name of the reaction, if set (otherwise returns the id)

    """
    r = model.select_one("listOfReactions").find(id=reaction_id)
    if "name" in r.attrs.keys() and r.attrs["name"]:
        return r.attrs["name"]
    else:
        return reaction_id
