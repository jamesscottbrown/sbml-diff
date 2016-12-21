from bs4 import BeautifulSoup
from accessor_functions import *
from generate_dot import *
from DiffObject import DiffObject
from rate_laws import *
from miriam import align_models
from tabulate import tabulate
import sys
import re


class SBMLDiff:

    def __init__(self, model_strings, model_names, generate_dot, align=False, cartoon=False, show_params=True, hide_rules=False):
        """

        Parameters
        ----------
        model_strings : a list, in which each element is an SBML model as a string
        model_names : names of each model (used as headings for the columns in table)
        generate_dot : instance of the GenerateDot class
        align : Boolean indicating whether to try to match using MIRIAM annotations as well as reaction/species id
        cartoon : Boolean indicating whether to draw transcription as a SBOLv promoter/CSD glyph, and hide degredation

        Returns
        -------
        models : list of models (each a bs4.BeautifulSoup object produced by parsing an SBML model)

        """

        self.model_strings = model_strings
        self.model_names = model_names
        self.generate_dot = generate_dot
        self.align = align
        self.cartoon = cartoon
        self.show_params = show_params
        self.hide_rules = hide_rules

        self.diff_object = DiffObject()

        self.models = map(lambda x: BeautifulSoup(x, 'xml'), self.model_strings)

        # Avoid need to search for reactions by id
        self.reactions = []
        for model in self.models:
            mr = {}
            reaction_list = model.select_one("listOfReactions")
            if reaction_list:
                for reaction in reaction_list.select("reaction"):
                    reaction_id = reaction.attrs["id"]
                    mr[reaction_id] = reaction
            self.reactions.append(mr)

        # avoid need to keep finding reactant compartments
        self.species_compartment = []
        self.initial_value = []
        for model in self.models:
            tmp = {}
            tmp_concentrations = {}

            species_list = model.select_one("listOfSpecies")
            if species_list:
                for species in species_list.select("species"):
                    species_id = species.attrs["id"]
                    compartment = species.attrs["compartment"]
                    tmp[id] = compartment

                    if "initialConcentration" in species.attrs:
                        tmp_concentrations[species_id] = species.attrs["initialConcentration"]

            self.species_compartment.append(tmp)
            self.initial_value.append(tmp_concentrations)

        # get initial parameter values
        self.initial_parameters = []
        for model_num, model in enumerate(self.models):
            for param in model.select("parameter"):
                if "id" not in param.attrs.keys():
                    continue
                param_id = param.attrs["id"]
                if "value" in param.attrs:
                    self.initial_value[model_num][param_id] = param.attrs["value"]

        # avoid need to search for reaction name
        self.reaction_name = []
        for model in self.models:
            tmp = {}

            reaction_list = model.select_one("listOfReactions")
            if reaction_list:
                for r in reaction_list.select("reaction"):
                    reaction_id = r.attrs["id"]
                    if "name" in r.attrs.keys() and r.attrs["name"]:
                        tmp[reaction_id] = r.attrs["name"]
                    else:
                        tmp[reaction_id] = reaction_id

            self.reaction_name.append(tmp)


        if self.cartoon:
            self.elided_list = []
            self.elided_reactions = []
            self.downstream_species = []
            self.find_downstream_species()

        self.modified_params = {}

    def check_model_supported(self):
        """
        Print an error message and quit if the file cannot be processed (because it contains user-defined functions, or is
        missing a list of species), rather than dumping a stack trace.
        """
        for model in self.models:

            if model.select_one('listOfReactions') and not model.select_one('listOfSpecies'):
                raise RuntimeError("Every model that includes a listOfReactions must include a listOfSpecies.")

            if not model.select_one('sbml') or 'xmlns' not in model.select_one('sbml').attrs.keys():
                raise RuntimeError("Every file must be an sbml model")

            if "level1" in model.select_one('sbml').attrs['xmlns']:
                raise RuntimeError("Every model must be in SBML level 2 or higher, since sbml-diff relies on id attributes")

    def print_rate_law_table(self, format="simple"):
        """
        Print a table of kineticLaws, in which rows correspond to reactions and columns to models.

        Parameters
        ----------
        format : a table format supported by tabulate (e.g. simple, html)
        """

        # get list of all reactions in all models
        reactions = []
        for model in self.models:
            reactions.extend(get_reactions(model))
        reactions = list(set(reactions))
        reactions.sort()

        rows = []
        for reaction_id in reactions:
            rates = [reaction_id]
            for model_num, model in enumerate(self.models):
                found_kinetic_law = False
                r = model.select_one("listOfReactions").find(id=reaction_id)
                if r:
                    kinetic_law = r.select_one("kineticLaw")
                    if kinetic_law:
                        math_tag = kinetic_law.select_one("math")
                        rates.append(convert_rate_law(math_tag))
                        found_kinetic_law = True

                if not found_kinetic_law:
                    rates.append("-")

                r = rates[1:]
                if r.count(r[0]) != len(r):
                    self.generate_dot.differences_found = True

            rows.append(rates)

        print tabulate(rows, ["Reaction"] + self.model_names, tablefmt=format)

    def compare_params(self, format="simple"):
        """
        Print a table of parameter values, in which rows correspond to reactions and columns to models.

        Parameters
        ----------
        format : a table format supported by tabulate (e.g. simple, html)
        """

        models = map(lambda x: BeautifulSoup(x, 'xml'), self.model_strings)

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

                p = row[1:]
                if p.count(p[0]) != len(p):
                    self.generate_dot.differences_found = True
            rows.append(row)

        print tabulate(rows, ["Parameter"] + self.model_names, tablefmt=format)

    def diff_events(self):
        """
        Compare all events between models.
        The id attribute is optional for event elements. For simplicity, we ignore ids even if they are present, so that
        two non-identical events between models are treated as entirely separate; it would be nicer if color of only
        those visual elements corresponding to what actually changed.
        """

        event_status = {}
        event_objects = {}
        has_id = {}

        for model_num, model in enumerate(self.models):
            event_list = model.select_one("listOfEvents")
            if not event_list:
                continue

            for event in event_list.select('event'):

                if 'id' not in event.attrs.keys():
                    event.attrs["id"] = str(hash(event))
                event_id = event.attrs["id"]

                if event_id not in event_status.keys():
                    event_status[event_id] = []

                event_status[event_id].append(model_num)
                event_objects[event_id] = event

        for event_id in event_objects.keys():
            self.diff_event_with_id(event_id, event_status[event_id])

    def diff_event_with_id(self, event_id, model_set):

        for model_num in model_set:
            species_ids = []
            species_list = self.models[model_num].select_one("listOfSpecies")
            if species_list:
                for s in species_list.select("species"):
                    species_ids.append(s.attrs["id"])

        # process trigger statement
        event_name = ""
        trigger_status = {}
        set_species_status = {}
        affects_value_status = {}
        modifier_arrows = {}
        trigger_expr = ""
        assignment_expr = {}
        trigger_param_status = {}
        assignment_param_arrows = {}

        for model_num in model_set:
            event = self.models[model_num].select_one('#' + event_id)

            # process model name
            if not event_name and "name" in event.attrs.keys():
                event_name = event.attrs["name"]

            # process trigger statements
            trigger = event.select_one("trigger")
            if trigger:
                for ci in trigger.select("ci"):
                    entity = ci.text.strip()
                    if entity in species_ids:
                        if entity not in trigger_status.keys():
                            trigger_status[entity] = []
                        trigger_status[entity].append(model_num)
                    else:
                        if entity not in trigger_param_status.keys():
                            trigger_param_status[entity] = []
                        trigger_param_status[entity].append(model_num)

                this_math_expr = trigger.select_one("math")
                if this_math_expr:
                    this_math_expr = convert_rate_law(this_math_expr)
                    if not trigger_expr or this_math_expr == trigger_expr:
                        trigger_expr = this_math_expr
                    else:
                        trigger_expr = "different"


            event_assignments = event.select("eventAssignment")
            if event_assignments:
                for event in event_assignments:
                    if isinstance(event, NavigableString):
                        continue

                    # arrow to species set
                    variable_id = event.attrs["variable"]
                    if variable_id in species_ids:

                        if variable_id not in set_species_status.keys():
                            set_species_status[variable_id] = []
                        set_species_status[variable_id].append(model_num)

                    elif self.show_params:
                        if variable_id not in set_species_status.keys():
                            set_species_status[variable_id] = []
                        set_species_status[variable_id].append(model_num)

                        if variable_id not in self.modified_params.keys():
                            self.modified_params[variable_id] = set()
                        self.modified_params[variable_id] = self.modified_params[variable_id].union(model_set)

                    # arrow from species affecting expression
                    math = event.select_one("math")
                    for ci in math.select("ci"):
                        species = ci.text.strip()
                        arrow_direction = categorise_interaction(math.parent, species, self.initial_value[model_num])
                        arrow = (species, arrow_direction, variable_id)

                        if species in species_ids:
                            if arrow not in modifier_arrows.keys():
                                modifier_arrows[arrow] = []

                            modifier_arrows[arrow].append(model_num)
                        else:
                            if arrow not in assignment_param_arrows.keys():
                                assignment_param_arrows[arrow] = []
                            assignment_param_arrows[arrow].append(model_num)

                    converted_math = convert_rate_law(math)
                    if variable_id not in assignment_expr.keys():
                        assignment_expr[variable_id] = converted_math
                    elif assignment_expr[variable_id] != converted_math:
                        assignment_expr[variable_id] = "different"

                    for arrow in modifier_arrows.keys():
                        if arrow not in affects_value_status.keys():
                            affects_value_status[arrow] = []
                        affects_value_status[arrow].append(model_num)

                    # arrow from param effe


        # record event node
        diff_event = self.diff_object.add_event()
        diff_event.set_event(event_id, event_name, model_set)

        for species in trigger_status:
            diff_event.add_trigger_species(species, event_id, trigger_status[species])
        diff_event.set_trigger(trigger_expr)

        for species in set_species_status:
            diff_event.add_set_species(species, event_id, assignment_expr[species], set_species_status[species])

        for arrow in affects_value_status:
            diff_event.add_event_affect_value_arrow(arrow[2], arrow[0], event_id, arrow[1], list(modifier_arrows[arrow]))

        for arrow in assignment_param_arrows:
            diff_event.add_assignment_param_arrow(arrow[2], arrow[0], event_id, arrow[1], list(assignment_param_arrows[arrow]))

        for param in trigger_param_status:
            diff_event.add_param(param, trigger_param_status[param], event_id)

    def diff_algebraic_rules(self):
        """
        Compare all algebraic rules between models.
        """

        rule_status = {}
        species_status = {}
        param_status = {}
        rate_laws = {}

        for model_num, model in enumerate(self.models):

            rule_list = model.select_one("listOfRules")
            if not rule_list:
                return

            # get details of each rule
            for rule in rule_list.select("algebraicRule"):

                # find species occurring in this rule
                species_ids = []
                species_list = model.select_one("listOfSpecies")
                species_in_rule = []
                params_in_rule = []
                if species_list:
                    for s in species_list.select("species"):
                        species_ids.append(s.attrs["id"])

                for ci in rule.select("ci"):
                    species_id = ci.string.strip()
                    if species_id in species_ids:
                        species_in_rule.append(species_id)
                    else:
                        params_in_rule.append(species_id)

                # Choose an id  to represent this rule
                if "metaid" in rule.attrs.keys():
                    rule_id = rule.attrs["metaid"]
                else:
                    rule_id = "assignmentRule" + "_".join(species_in_rule)

                # record that model contained this rule
                if rule_id not in rule_status.keys():
                    rule_status[rule_id] = []
                    species_status[rule_id] = {}
                    param_status[rule_id] = {}
                rule_status[rule_id].append(model_num)

                # record species
                for species_id in species_in_rule:
                    if species_id not in species_status[rule_id].keys():
                        species_status[rule_id][species_id] = []

                    species_status[rule_id][species_id].append(model_num)

                # record params
                for param_id in params_in_rule:
                    if param_id not in param_status[rule_id].keys():
                        param_status[rule_id][species_id] = []
                    param_status[rule_id][species_id].append(model_num)

                # record math expression
                rate_law = rule.select_one("math")

                if not rate_law:
                    continue

                if rate_law and not rule_id in rate_laws.keys():
                    rate_laws[rule_id] = rate_law
                if rule_id in rate_laws.keys() and rate_law and rate_laws[rule_id] != rate_law:
                    rate_laws[rule_id] = "different"

        # produce output
        for rule_id in rule_status.keys():
            rule = self.diff_object.add_rule()

            for species_id in species_status[rule_id]:
                rule.add_algebraic_arrow(species_status[rule_id][species_id], rule_id, species_id)

            for param_id in param_status[rule_id]:
                rule.add_parameter_rule(param_status[rule_id][param_id], rule_id, param_id, 'none')

            converted_rate_law = ""
            if rule_id in rate_laws.keys() and rate_laws[rule_id] != "different":
                converted_rate_law = convert_rate_law(rate_laws[rule_id])

            rule.set_rule(rule_status[rule_id], rule_id, rate_laws, converted_rate_law)

    def diff_rules(self):
        """
        Compare all (rate or assignment) rules between models.
        """
        rule_status = {}
        for model_num, model in enumerate(self.models):
            rule_targets = get_variables_set_by_rules(model)

            for rule_target in rule_targets:
                species_list = model.select_one('listOfSpecies')
                if not species_list or not species_list.find(id=rule_target):
                    if rule_target not in self.modified_params.keys():
                        self.modified_params[rule_target] = set()
                    self.modified_params[rule_target].add(model_num)

                if rule_target not in rule_status.keys():
                    rule_status[rule_target] = set()

                rule_status[rule_target].add(model_num)

        for rule_target in rule_status:
            self.diff_rule(rule_target)

    def diff_rule(self, target_id):
        """
        Compare a single rule between models.

        Parameters
        ----------
        target_id : id of the species affected by this rule
        """
        # if a reaction is shared, we need to consider whether its products, reactants and rate law are also shared

        # 'modifiers' appear in the math expression of a rule that sets 'target'
        # a rule has only one target, whereas reaction may have multiple products

        modifier_status = {}
        parameter_status = {}
        target_status = set()
        rate_laws = ""

        rule = self.diff_object.add_rule()
        for model_num, model in enumerate(self.models):
            modifiers, compartment, rate_law = get_rule_details(model, target_id, self.species_compartment[model_num])

            if not rate_law:
                continue

            if rate_law and not rate_laws:
                rate_laws = rate_law
            if rate_laws and rate_law and rate_laws != rate_law:
                rate_laws = "different"

            for modifier in modifiers:
                arrow_direction = categorise_interaction(rate_law.parent, modifier, self.initial_value[model_num])
                arrow = (modifier, arrow_direction)
                if arrow not in modifier_status.keys():
                    modifier_status[arrow] = set()
                modifier_status[arrow].add(model_num)

            entities = rate_law.select("ci")
            for entity in entities:
                param = entity.string.strip()

                # check a param rather than species
                if param in self.species_compartment[model_num].keys():
                    continue

                arrow_direction = categorise_interaction(rate_law.parent, param, self.initial_value[model_num])
                arrow = (param, arrow_direction)

                if arrow not in parameter_status.keys():
                    parameter_status[arrow] = set()
                parameter_status[arrow].add(model_num)

            # for target in targets:
            target_status.add(model_num)

        # modifier arrows
        for arrow in modifier_status:
            model_set = list(modifier_status[arrow])
            rule.add_modifier_arrow(model_set, target_id, arrow[0], arrow[1])

        # parameter arrows
        for arrow in parameter_status:
            model_set = list(parameter_status[arrow])
            rule.add_parameter_rule(model_set, target_id, arrow[0], arrow[1])

        # target arrows
        model_set = list(target_status)
        species_list = self.models[model_set[0]].select_one('listOfSpecies')
        if self.show_params or (species_list and species_list.find(id=target_id)):
            rule.add_target_arrow(model_set, target_id)

        # rate law
        converted_rate_law = ""
        if rate_laws and rate_laws != "different":
            converted_rate_law = convert_rate_law(rate_laws)

        rule.set_rule(model_set, target_id, rate_laws, converted_rate_law)

    def diff_reactions(self):
        """
        Compare all reactions between models.
        """

        reaction_status = {}
        for model_num, model in enumerate(self.models):
            reactions = get_reactions(model)

            for reaction in reactions:
                if reaction not in reaction_status.keys():
                    reaction_status[reaction] = set()

                reaction_status[reaction].add(model_num)

        for reaction_id in reaction_status:
            self.diff_reaction(reaction_id)

    def diff_reaction(self, reaction_id):
        """
        Compare a single reaction between models.

        Parameters
        ----------
        reaction_id : id of the reaction
        """

        # We need to consider whether the reaction's products, reactants and rate law are shared
        reaction_model_set = set()
        fast_model_set = set()
        irreversible_model_set = set()

        reactant_status = {}
        product_status = {}
        parameter_status = {}

        compartment = ""
        rate_laws = ""
        rate_law_found = False

        reactant_stoichiometries = {}
        product_stoichiometries = {}

        transcription_reaction = False
        ever_drawn = False

        for model_num, model in enumerate(self.models):
            if reaction_id not in self.reactions[model_num].keys():
                continue
            reaction = self.reactions[model_num][reaction_id]

            reactants, products, compartment, rate_law, rs, ps = get_reaction_details(model, reaction, self.species_compartment[model_num])

            # Skip processing reaction if it should not be drawn for this model
            show_reaction = True
            if self.cartoon:
                if reaction in self.elided_reactions[model_num]:
                    show_reaction = False

                # If a reaction has only one reaction, and it is an elided species (e.g. translation, mRNA degredation), do not print anything
                if len(reactants) == 1 and reactants[0] in self.elided_list[model_num]:
                    show_reaction = False
                # TODO: what if only a modifier

                # Hide all degredaion
                if len(reactants) == 1 and len(products) == 0:
                    show_reaction = False

            if not show_reaction:
                continue

            if self.cartoon and "sboTerm" in reaction.attrs.keys() and \
                    reaction.attrs['sboTerm'] in ["SBO:0000183", "SBO:0000589"]:
                transcription_reaction = True

            # only perform comparison between models in which this reaction actually occurs
            if not reactants and not products and not compartment and not rate_law and not rs and not ps:
                continue

            reaction_model_set.add(model_num)

            if "fast" in reaction.attrs.keys() and reaction.attrs["fast"] in ['1', 'true']:
                fast_model_set.add(model_num)
            if "reversible" in reaction.attrs.keys() and reaction.attrs["reversible"] in ['0', 'false']:
                irreversible_model_set.add(model_num)

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

            if not rate_law_found:
                rate_laws = rate_law
                rate_law_found = True
            if rate_law_found and re.sub('\s', '', str(rate_law)) != re.sub('\s', '', str(rate_laws)):
                rate_laws = "different"

            for reactant in reactants:
                if reactant not in reactant_status.keys():
                    reactant_status[reactant] = set()
                reactant_status[reactant].add(model_num)

            if rate_law:
                entities = rate_law.select("ci")
                for entity in entities:
                    param = entity.string.strip()

                    # check a param rather than species
                    if param in self.species_compartment[model_num].keys():
                        continue

                    arrow_direction = categorise_interaction(rate_law.parent, param, self.initial_value[model_num])
                    arrow = (param, arrow_direction)

                    if arrow not in parameter_status.keys():
                        parameter_status[arrow] = set()
                    parameter_status[arrow].add(model_num)

            for product in products:
                # if producing something that's been elided, adjust arrows to point ot downstream species
                if self.cartoon and product in self.elided_list[model_num]:
                    new_product = self.downstream_species[model_num][product]
                    product_stoichiometries[new_product] = product_stoichiometries[product]
                    product = new_product

                if product not in product_status.keys():
                    product_status[product] = set()
                product_status[product].add(model_num)

            parent_model = model_num
            ever_drawn = True

        # If reaction should not be drawn for any models, return now
        if not ever_drawn:
            return ""

        reaction = self.diff_object.add_reaction(compartment)

        # reactant arrows
        for reactant_num, reactant in enumerate(reactant_status):
            model_set = list(reactant_status[reactant])
            reaction.add_reactant_arrow(model_set, reaction_id, reactant, reactant_stoichiometries[reactant])

        # product arrows
        for product_num, product in enumerate(product_status):
            model_set = list(product_status[product])
            if transcription_reaction:
                reaction.add_transcription_product_arrow(model_set, reaction_id, product, product_stoichiometries[product])
            else:
                reaction.add_product_arrow(model_set, reaction_id, product, product_stoichiometries[product])

        # parameter arrows
        for arrow in parameter_status:
            model_set = list(parameter_status[arrow])
            reaction.add_parameter_arrow(model_set, reaction_id, arrow[0], arrow[1])

        # rate law
        reaction_name = self.reaction_name[parent_model][reaction_id]

        converted_rate_law = ""
        if rate_laws and rate_laws != "different":
            converted_rate_law = convert_rate_law(rate_laws)

        reaction.set_reaction(reaction_model_set, reaction_id, rate_laws, reaction_name, converted_rate_law,
                                     fast_model_set, irreversible_model_set, product_status, transcription_reaction)

    def find_downstream_species(self):
        """
        Identifies reactions which should be elided in cartoon mode. A reaction should be elided if:
        - it has sboTerm SBO:0000184 (translation)
        - it has exactly one modifier or reactant species, and this species does not feature as a reactant or modifier
        species of any reaction that is not translation or degredation
        - it has exactly one product
        - it does not appear in more than one model, unless it has the same kineticLaw in each (as eliding the reaction
        would hide this difference)
        """

        for model_num, model in enumerate(self.models):
            # Only elide reactions with sboTerm corresponding to translation, and only one reactant/modifier species

            self.elided_list.append([])
            self.elided_reactions.append([])
            self.downstream_species.append({})


            # first, form a list of species that cannot safely be elided, because they are a reactant or modifier in a
            # reaction other than degredation or translation
            non_intermediates = []
            for reaction in model.select('reaction'):

                # skip degredation or translation reactions
                if "sboTerm" in reaction.attrs.keys() and reaction.attrs["sboTerm"] in ["SBO:0000184", "SBO:0000179"]:
                    continue

                reactant_list = reaction.select_one("listOfReactants")
                if reactant_list:
                    for reactant in reactant_list.select("speciesReference"):
                        non_intermediates.append(reactant["id"])

                modifier_list = reaction.select_one("listOfModifiers")
                if modifier_list:
                    for r in modifier_list.select("modifierSpeciesReference"):
                        non_intermediates.append(r["species"])

            # Now loop through reactions, identifying those that should be elided
            for reaction in model.select('reaction'):

                if "sboTerm" not in reaction.attrs.keys() or reaction.attrs["sboTerm"] != "SBO:0000184":
                    continue

                # if reaction has different kineticLaw in different models, don't elide it
                rate_laws = ""
                for m in self.models:
                    r = m.select_one("listOfReactions").find(id=reaction["id"])
                    if not r:
                        continue

                    rate_law = r.select_one("kineticLaw").select_one("math")
                    if rate_law and not rate_laws:
                        rate_laws = rate_law
                    elif rate_laws and rate_law and rate_laws != rate_law:
                        rate_laws = "different"
                        break

                if rate_laws == "different":
                    continue

                reactants_and_modifier_species = []

                # Check exactly one modifier/reactant
                modifier_list = reaction.select_one("listOfModifiers")
                if modifier_list:
                    for r in modifier_list.select("modifierSpeciesReference"):
                        reactants_and_modifier_species.append(r["species"])

                reactant_list = reaction.select_one("listOfReactants")
                if reactant_list:
                    for r in reactant_list.select("speciesReference"):
                        reactants_and_modifier_species.append(r["species"])

                if len(reactants_and_modifier_species) != 1:
                    continue

                species_to_elide = reactants_and_modifier_species[0]

                if species_to_elide in non_intermediates:
                    continue

                # check exactly one product (other than reactant, in case reaction is modelled as mRNA -> mRNA + protein)
                product_species = []
                product_list = reaction.select_one("listOfProducts")
                if product_list:
                    for p in product_list.select("speciesReference"):
                        product_id = p["species"]
                        if product_id != species_to_elide:
                            product_species.append(product_id)

                if len(product_species) != 1:
                    continue

                self.elided_list[model_num].append(species_to_elide)
                self.elided_reactions[model_num].append(reaction)
                self.downstream_species[model_num][species_to_elide] = product_species[0]

    def diff_compartment(self, compartment_id):
        """
        Print DOT output comparing a single compartment between models

        Parameters
        ----------
        compartment_id : the id of a compartment
        reaction_strings : list of strings, each a DOT statement describing the nodes representing reaction
        rule_strings : list of strings, each a DOT statement describing the nodes representing assignmentRules/rateRules
        algebraic_rule_strings : list of strings, each a DOT statement describing the nodes representing algebraicRules
        """

        # For each species, find set of models containing it
        species_status = {}
        is_boundary_species = {}
        for model_num, model in enumerate(self.models):
            for species in get_species(model, compartment_id):

                if species not in species_status.keys():
                    species_status[species] = set()
                species_status[species].add(model_num)

                s = model.select_one("listOfSpecies").find(id=species)
                is_boundary = ""
                if "boundaryCondition" in s.attrs.keys():
                    is_boundary = s.attrs["boundaryCondition"]

                if species not in is_boundary_species.keys():
                    is_boundary_species[species] = is_boundary
                elif is_boundary_species[species] != is_boundary:
                    is_boundary_species[species] = '?'

        for species in species_status:
            parent_model_index = list(species_status[species])[0]
            parent_model = self.models[parent_model_index]
            species_name = get_species_name(parent_model, species)

            # In cartoon mode, don't draw species node unless there is >=1 model in which it is present and not elided
            draw = True
            if self.cartoon:
                draw = False
                for model_num in species_status[species]:
                    if species not in self.elided_list[model_num]:
                        draw = True
            if draw:
                self.diff_object.add_species(compartment_id, species_status[species], is_boundary_species[species], species, species_name)

        # for each regulatory interaction - (reactant, reaction, effect direction) tuple - find set of models containing it
        arrow_status = {}
        for model_num, model in enumerate(self.models):
            if self.cartoon:
                arrows = get_regulatory_arrow(model, compartment_id, self.reactions[model_num], self.species_compartment[model_num], self.initial_value[model_num], elided_reactions=self.elided_reactions[model_num])
            else:
                arrows = get_regulatory_arrow(model, compartment_id, self.reactions[model_num], self.species_compartment[model_num], self.initial_value[model_num])

            for ind, arrow in enumerate(arrows):
                if arrow not in arrow_status.keys():
                    arrow_status[arrow] = set()
                arrow_status[arrow].add(model_num)

        for ind, arrow in enumerate(arrow_status):
            self.diff_object.add_regulatory_arrow(compartment_id, arrow_status[arrow], arrow[0], arrow[1], arrow[2])

    def diff_models(self):
        """
        Print DOT output comparing SBML models
        """

        self.check_model_supported()
        self.models = map(lambda x: inline_all_functions(x), self.models)

        if self.align:
            align_models(self.models)

        self.diff_reactions()

        if not self.hide_rules:
            self.diff_rules()
            self.diff_algebraic_rules()


        compartment_ids = set()
        for model_num, model in enumerate(self.models):
            for compartment in model.select('compartment'):
                if "id" in compartment.attrs.keys():
                    compartment_ids.add(compartment.attrs["id"])

        for compartment_id in compartment_ids:
            self.diff_compartment(compartment_id)

        self.diff_events()
        if self.show_params:
            self.draw_modified_params()

        # actually print the results of comparison
        self.generate_dot.generate_dot(self.diff_object)

    def abstract_model(self, model, model_num):
        """
        For each pair of species in a model, determine if there is an interaction, and if so classify it.
        If a species is a reactant, it is not considered to interact with itself through that reaction, since any species
        increases the rate of its own degredation

        Parameters
        ----------
        model : bs4.BeautifulSoup object produced by parsing an SBML model

        model_num : index of the model being abstracted

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
        for reaction_id in reactions:
            reaction = self.reactions[model_num][reaction_id]
            reactant_list, product_list, compartment, rate_law, _, _ = get_reaction_details(model, reaction, self.species_compartment[model_num])

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

                    effect = categorise_interaction(rate_law.parent, modifier, self.initial_value[model_num])
                    if effect == "monotonic_increasing":
                        interactions[modifier][reactant].add("increase-degredation")
                    elif effect == "monotonic_decreasing":
                        interactions[modifier][reactant].add("decrease-degredation")

                for product in product_list:
                    effect = categorise_interaction(rate_law.parent, modifier, self.initial_value[model_num])
                    if effect == "monotonic_increasing":
                        interactions[modifier][product].add("increase-production")
                    elif effect == "monotonic_decreasing":
                        interactions[modifier][product].add("decrease-production")

        return interactions, species

    # TODO: compartments!
    def diff_abstract_models(self, ignored_species, elided_species):
        """
        Print DOT output comparing SBML models after abstraction (by abstract_model(model))

        Parameters
        ----------
        ignored_species : list of species to be simply removed
        elided_species : list of species to be removed, with interactions targeting them appropriately moved downstream

        """
        if not ignored_species:
            ignored_species = []
        if not elided_species:
            elided_species = []

        if self.align:
            align_models(self.models)

        effect_types = ["increase-degredation", "decrease-degredation", "increase-production", "decrease-production"]

        # Construct abstracted version of each model
        abstracted_model = [] * len(self.models)
        species_list = set()
        models_containing_species = {}
        is_boundary_species = {}

        for model_num, model in enumerate(self.models):
            abstract, species = self.abstract_model(model, model_num)

            abstracted_model.append(abstract)
            species_list = species_list.union(species)

            for s in species:
                if s not in models_containing_species.keys():
                    models_containing_species[s] = set()
                models_containing_species[s].add(model_num)

                species_object = model.select_one("listOfSpecies").find(id=s)
                is_boundary = ""
                if "boundaryCondition" in species_object.attrs.keys():
                    is_boundary = species_object.attrs["boundaryCondition"]

                if s not in is_boundary_species.keys():
                    is_boundary_species[s] = is_boundary
                elif is_boundary_species[s] != is_boundary:
                    is_boundary_species[s] = '?'

        species_list = species_list.difference(ignored_species)
        retained_species = species_list.difference(elided_species)

        self.generate_dot.print_header()

        for s in retained_species:
            model_num = list(models_containing_species[s])[0]
            species_name = get_species_name(self.models[model_num], s)
            self.generate_dot.print_species_node(models_containing_species[s], is_boundary_species[s], s, species_name)

        # Construct interactions[modifier][species][type] = set of model_numbers
        interactions = {}
        for s1 in species_list:
            interactions[s1] = {}
            for s2 in species_list:
                interactions[s1][s2] = {}
                for effect in effect_types:
                    interactions[s1][s2][effect] = set()

        for model_num, model in enumerate(self.models):
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
            interactions = self.elide(species_list, effect_types, interactions, elided_species)

        for modifier in retained_species:
            for species in retained_species:
                for effect_type in effect_types:
                    model_list = interactions[modifier][species][effect_type]
                    self.generate_dot.print_abstracted_arrow(model_list, modifier, species, effect_type)

        self.generate_dot.print_footer()

    def elide(self, species_list, effect_types, interactions, elided_species):
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
        effect_types : list of the possible effect types

        Returns
        -------
         a modified interactions structure
        """
        elided_species = set(elided_species).intersection(species_list)
        for model_num, model in enumerate(self.models):

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

    def draw_modified_params(self):
        for param_id in self.modified_params.keys():
            model_set = list(self.modified_params[param_id])
            name = param_id
            param = self.models[model_set[0]].find(id=param_id)
            if "name" in param.attrs.keys():
                name = param.attrs["name"]
            self.diff_object.add_param_node(param_id, name, model_set)
