from bs4 import BeautifulSoup
from accessor_functions import *
from generate_dot import *
from rate_laws import *
from tabulate import tabulate


def print_rate_law_table(model_strings, model_names):
    """
    Print a table of kineticLaws, in which rows correspond to reactions and columns to models.

    Parameters
    ----------
    model_strings : a list, in which each element is an SBML model as a string
    model_names : headings for the columns
    """

    # get list of all reactions in all models
    reactions = []
    models = map(lambda x: BeautifulSoup(x, 'xml'), model_strings)
    for model in models:
        reactions.extend(get_reactions(model))
    reactions = list(set(reactions))
    reactions.sort()

    rows = []
    for reaction_id in reactions:
        rates = [reaction_id]
        for model_num, model in enumerate(models):
            r = model.select_one("listOfReactions").find(id=reaction_id)
            if r:
                math_tag = r.select_one("kineticLaw").select_one("math")
                rates.append(convert_rate_law(math_tag))
            else:
                rates.append("-")
        rows.append(rates)

    print tabulate(rows, ["Reaction"] + model_names)


def compare_params(model_strings, model_names):
    """
    Print a table of parameter values, in which rows correspond to reactions and columns to models.

    Parameters
    ----------
    model_strings : a list, in which each element is an SBML model as a string
    model_names : headings for the columns
    """

    models = map(lambda x: BeautifulSoup(x, 'xml'), model_strings)

    param_value = {}
    for model_num, model in enumerate(models):
        param_ids, param_values = get_params(model)

        for param_id in param_ids:

            if param_id not in param_value.keys():
                param_value[param_id] = {}
            param_value[param_id][model_num] = param_values[param_id]

    rows = []
    for param_id in param_value.keys():
        row = [param_id]
        for model_num, model in enumerate(models):
            if model_num in param_value[param_id].keys():
                row.append(param_value[param_id][model_num])
            else:
                row.append("-")
        rows.append(row)

    print tabulate(rows, ["Parameter"] + model_names)


def diff_rules(models, generate_dot):
    """
    Compare all rules between models. Returns a list of the DOT statements to draw each rule node, but directly
    prints the DOT statements for the corresponding arrows.

    Parameters
    ----------
    models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)
    generate_dot : instance of the GenerateDot class

    Returns
    -------
    a list of the DOT statements to draw each rule node

    """
    rule_status = {}
    for model_num, model in enumerate(models):
        rule_targets = get_species_set_by_rules(model)

        for rule_target in rule_targets:
            if rule_target not in rule_status.keys():
                rule_status[rule_target] = set()

            rule_status[rule_target].add(model_num)

    rule_strings = {}

    for rule_target in rule_status:
        model_set = list(rule_status[rule_target])
        inputs, compartment, rate_law = get_rule_details(models[model_set[0]], rule_target)

        if compartment not in rule_strings.keys():
            rule_strings[compartment] = []

        rule_string = diff_rule(models, rule_target, generate_dot)
        rule_strings[compartment].append(rule_string)

    return rule_strings


def diff_rule(models, target_id, generate_dot):
    """
    Compare a single rule between models. Returns the DOT statement to draw the node for the rule, but directly
    prints the DOT statements for the corresponding arrows. This is to allow the rule node to be drawn in the correct
    compartment.

    Parameters
    ----------
    models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)
    target_id : id of the species affected by this rule
    generate_dot : instance of the GenerateDot class

    Returns
    -------
    a DOT statement describing the nodes representing this rule
    """
    # if a reaction is shared, we need to consider whether its products, reactants and rate law are also shared

    # 'modifiers' appear in the math expression of a rule that sets 'target'
    # a rule has only one target, whereas reaction may have multiple products

    modifier_status = {}
    target_status = set()
    rate_laws = ""

    for model_num, model in enumerate(models):
        modifiers, compartment, rate_law = get_rule_details(model, target_id)

        if not rate_law:
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
        generate_dot.print_rule_modifier_arrow(model_set, target_id, modifier)

    # target arrows
    model_set = list(target_status)
    generate_dot.print_rule_target_arrow(model_set, target_id)

    # rate law
    converted_rate_law = ""
    if rate_laws and rate_laws != "different":
        converted_rate_law = convert_rate_law(rate_laws)

    return generate_dot.print_rule_node(model_set, target_id, rate_laws, converted_rate_law)


def diff_reactions(models, generate_dot):
    """
    Compare all reactions between models. Returns a list of the DOT statements to draw each reaction node, but directly
    prints the DOT statements for the corresponding arrows.

    Parameters
    ----------
    models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)
    generate_dot : instance of the GenerateDot class

    Returns
    -------
     a list of the DOT statements to draw each reaction node
    """

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
        reactant_list, product_list, compartment, rate_law, _, _ = get_reaction_details(models[model_set[0]], reaction_id)

        if compartment not in reaction_strings.keys():
            reaction_strings[compartment] = []

        reaction_string = diff_reaction(models, reaction_id, generate_dot)
        reaction_strings[compartment].append(reaction_string)

    return reaction_strings


def diff_reaction(models, reaction_id, generate_dot):
    """
    Compare a single reaction between models. Returns the DOT statement to draw the node for the reaction, but directly
    prints the DOT statements for the corresponding arrows. This is to allow the reaction node to be drawn in the correct
    compartment.

    Parameters
    ----------
    models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)
    reaction_id : id of the reaction
    generate_dot : instance of the GenerateDot class

    Returns
    -------
     a DOT statement describing the node representing this reaction
    """

    # We need to consider whether the reaction's products, reactants and rate law are shared

    reaction_model_set = set()
    reactant_status = {}
    product_status = {}
    rate_laws = ""

    reactant_stoichiometries = {}
    product_stoichiometries = {}

    for model_num, model in enumerate(models):
        reactants, products, compartment, rate_law, rs, ps = get_reaction_details(model, reaction_id)

        # only perform comparison between models in which this reaction actually occurs
        if not reactants and not products and not compartment and not rate_law and not rs and not ps:
            continue

        reaction_model_set.add(model_num)

        # check if any stoichiometry values change between models
        # record the stoichiometry of every reactant or product associated with this reaction in any model
        # if a reactant/product has 2 or more storichiometries between the models, record it as '?'
        for ind, stoich in enumerate(rs):
            if reactants[ind] not in reactant_stoichiometries.keys():
                reactant_stoichiometries[reactants[ind]] = stoich
            elif stoich != reactant_stoichiometries[reactants[ind]]:
                reactant_stoichiometries[reactants[ind]] = '?'

        for ind, stoich in enumerate(ps):
            if products[ind] not in product_stoichiometries.keys():
                product_stoichiometries[products[ind]] = stoich
            elif stoich != product_stoichiometries[products[ind]]:
                product_stoichiometries[products[ind]] = '?'

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
    for reactant_num, reactant in enumerate(reactant_status):
        model_set = list(reactant_status[reactant])
        generate_dot.print_reactant_arrow(model_set, reaction_id, reactant, reactant_stoichiometries[reactant])

    # product arrows
    for product_num, product in enumerate(product_status):
        model_set = list(product_status[product])
        generate_dot.print_product_arrow(model_set, reaction_id, product, product_stoichiometries[product])

    # rate law
    parent_model = models[model_set[0]]
    reaction_name = get_reaction_name(parent_model, reaction_id)

    converted_rate_law = ""
    if rate_laws and rate_laws != "different":
        converted_rate_law = convert_rate_law(rate_laws)

    return generate_dot.print_reaction_node(reaction_model_set, reaction_id, rate_laws, reaction_name, converted_rate_law)


def diff_compartment(compartment_id, models, reaction_strings, rule_strings, generate_dot):
    """
    Print DOT output comparing a single compartment between models

    Parameters
    ----------
    compartment_id : the id of a compartment
    models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)
    reaction_strings : list of strings, each a DOT statement describing the nodes representing reaction
    rule_strings : list of strings, each a DOT statement describing the nodes representing rules
    generate_dot : instance of the GenerateDot class
    """
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


def diff_models(model_strings, generate_dot):
    """
    Print DOT output comparing SBML models

    Parameters
    ----------
    model_strings : a list, in which each element is an SBML model as a string

    generate_dot : instance of the GenerateDot class

    """
    models = map(lambda x: BeautifulSoup(x, 'xml'), model_strings)

    generate_dot.print_header()

    reaction_strings = diff_reactions(models, generate_dot)
    rule_strings = diff_rules(models, generate_dot)

    if "NONE" in reaction_strings.keys():
        print "\n".join(reaction_strings["NONE"])
    if "NONE" in rule_strings.keys():
        print "\n".join(rule_strings["NONE"])

    compartment_ids = set()
    for model_num, model in enumerate(models):
        for compartment in model.select('compartment'):
            compartment_ids.add(compartment.attrs["id"])

    for compartment_id in compartment_ids:
        diff_compartment(compartment_id, models, reaction_strings, rule_strings, generate_dot)

    generate_dot.print_footer()


def abstract_model(model):
    """
    For each pair of species in a model, determine if there is an interaction, and if so classify it.
    If a species is a reactant, it is not considered to interact with itself through that reaction, since any species
    increases the rate of its own degredation

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model

    Returns
    -------
    interactions : a 2D list, in which entries interactions[modifier][product] indicate the effect of the species with
        id modifier on the species with id product - one of "increase-degredation", "decrease-degredation",
        "increase-production", or "decrease-production"

    species : id of each species in the model
    """
    # Get list of species
    species = set()
    for compartment in model.select('compartment'):
        compartment_id = compartment.attrs["id"]
        species = species.union(get_species(model, compartment_id))

    interactions = {}
    for modifier in species:
        interactions[modifier] = {}
        for target in species:
            interactions[modifier][target] = set()

    reactions = get_reactions(model)
    for reaction in reactions:
        reactant_list, product_list, compartment, rate_law, _, _ = get_reaction_details(model, reaction)

        # Identify all species that appear in kineticLaw
        modifiers = []
        for ci in rate_law.findAll("ci"):
            name = ci.text.strip()
            if name in species:
                modifiers.append(name)
        modifiers = set(modifiers)

        for modifier in modifiers:
            for reactant in reactant_list:

                # Any species increases the rate of its own degredation, so ignore this
                if reactant == modifier:
                    continue

                effect = categorise_interaction(rate_law.parent, modifier)
                if effect == "monotonic_increasing":
                    interactions[modifier][reactant].add("increase-degredation")
                elif effect == "monotonic_decreasing":
                    interactions[modifier][reactant].add("decrease-degredation")

            for product in product_list:
                effect = categorise_interaction(rate_law.parent, modifier)
                if effect == "monotonic_increasing":
                    interactions[modifier][product].add("increase-production")
                elif effect == "monotonic_decreasing":
                    interactions[modifier][product].add("decrease-production")

    return interactions, species


# TODO: compartments!
def diff_abstract_models(model_strings, generate_dot, ignored_species=False, elided_species=False):
    """
    Print DOT output comparing SBML models after abstraction (by abstract_model(model))

    Parameters
    ----------
    model_strings : a list, in which each element is an SBML model as a string
    generate_dot : instance of the GenerateDot class
    ignored_species : list of species to be simply removed
    elided_species : list of species to be removed, with interactions targeting them appropriately moved downstream
    """

    if not ignored_species:
        ignored_species = []
    if not elided_species:
        elided_species = []

    models = map(lambda x: BeautifulSoup(x, 'xml'), model_strings)

    effect_types = ["increase-degredation", "decrease-degredation", "increase-production", "decrease-production"]

    # Construct abstracted version of each model
    abstracted_model = [] * len(models)
    species_list = set()
    models_containing_species = {}

    for model_num, model in enumerate(models):
        abstract, species = abstract_model(model)

        abstracted_model.append(abstract)
        species_list = species_list.union(species)

        for s in species:
            if s not in models_containing_species.keys():
                models_containing_species[s] = set()
            models_containing_species[s].add(model_num)

    species_list = species_list.difference(ignored_species)
    retained_species = species_list.difference(elided_species)

    print "digraph comparison {"

    for s in retained_species:
        model_num = list(models_containing_species[s])[0]
        species_name = get_species_name(models[model_num], s)
        generate_dot.print_species_node(models_containing_species[s], s, species_name)

    # Construct interactions[modifier][species][type] = set of model_numbers
    interactions = {}
    for s1 in species_list:
        interactions[s1] = {}
        for s2 in species_list:
            interactions[s1][s2] = {}
            for effect in effect_types:
                interactions[s1][s2][effect] = set()

    for model_num, model in enumerate(models):
        for modifier in species_list:
            if model_num not in models_containing_species[modifier]:
                continue

            for species in species_list:
                if model_num not in models_containing_species[species]:
                    continue

                effects = abstracted_model[model_num][modifier][species]
                for effect_type in effects:
                    interactions[modifier][species][effect_type].add(model_num)

    if elided_species:
        interactions = elide(interactions, elided_species, species_list, models, effect_types)

    for modifier in retained_species:
        for species in retained_species:
            for effect_type in effect_types:
                model_list = interactions[modifier][species][effect_type]
                generate_dot.print_abstracted_arrow(model_list, modifier, species, effect_type)

    print "}"


def elide(interactions, elided_species, species_list, models, effect_types):
    """
    Removes certain species from a model, transfering interactions that target them onto the species that they produce.
    Intended for case of an intermediate that causes production of a downstream species (e.g. mRNA, which causes
    production of a protein)

    Parameters
    ----------
    interactions : interactions[modifier][target][effect_type] is set of model numbers for which species with id
                    modifier has effect effect_type on species with id target
    elided_species : list containing the ide of each species to elide
    species_list : list containing id of every species
    models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)
    effect_types : list of the possible effect types

    Returns
    -------
     a modified interactions structure
    """
    elided_species = set(elided_species).intersection(species_list)
    for model_num, model in enumerate(models):

        # For each elided species,
        for s in elided_species:

            # find the 'downstream' species (eg. the protein produced from mRNA)
            downstream = False
            for s2 in species_list:
                if model_num in interactions[s][s2]["increase-production"]:
                    downstream = s2

            if not downstream:
                continue

            # Transfer interactions targeting the elided species to the downstream species
            for regulator in species_list:
                for effect_type in effect_types:
                    if model_num in interactions[regulator][s][effect_type]:
                        interactions[regulator][downstream][effect_type].add(model_num)

    # Then remove the elided species
    for s in elided_species:
        interactions.pop(s)

        for s2 in species_list:
            if s2 in interactions.keys():
                interactions[s2].pop(s)

    return interactions
