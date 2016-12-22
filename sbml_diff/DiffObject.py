import collections


class DiffObject:
    def __init__(self):
        self.compartments = {}
        self.add_compartment("NONE")
        self.events = []  # should probably be moved into compartment ?
        self.param_nodes = []

    def add_compartment(self, compartment_id):
        self.compartments[compartment_id] = DiffCompartment()
        return self.compartments[compartment_id]

    def check_compartment_exists(self, compartment):
        if compartment not in self.compartments.keys():
            self.add_compartment(compartment)
        return self.compartments[compartment]

    def add_event(self):
        new_event = DiffEvent()
        self.events.append(new_event)
        return new_event

    def add_param_node(self, variable_id, variable_name, model_set):
        self.param_nodes.append({"variable_id": variable_id, "variable_name": variable_name, "model_set": model_set})


class DiffCompartment:
    def __init__(self):
        self.species = {}
        self.regulatory_arrows = DiffElement()
        self.reactions = []
        self.rules = []

    def add_species(self, species_id, is_boundary, species_name, elided, model_num):

        if species_id not in self.species.keys():
            self.species[species_id] = DiffElement()

        self.species[species_id].add({"species_id": species_id, "is_boundary": is_boundary, "species_name": species_name,
                          "elided": elided}, model_num)

    def add_regulatory_arrow(self, arrow_source, arrow_target, arrow_direction, model_num):
        self.regulatory_arrows.add({"arrow_source": arrow_source, "arrow_target": arrow_target,
                                    "arrow_direction": arrow_direction}, model_num)

    def add_reaction(self):
        new_reaction = DiffReaction()
        self.reactions.append(new_reaction)
        return new_reaction

    def add_rule(self, rule_id):
        new_rule = DiffRule(rule_id)
        self.rules.append(new_rule)
        return new_rule


class DiffEventAssignment:
    def __init__(self):
        self.affect_value_arrows = DiffElement()
        self.affect_value_param_arrows = DiffElement()
        self.math_expr = DiffElement()


class DiffEvent:
    def __init__(self):
        self.event = {}
        self.trigger_arrows = DiffElement()
        self.assignments = {}
        self.trigger_math = DiffElement()
        self.trigger_params = DiffElement()

    def check_target_exists(self, target):
        if target not in self.assignments.keys():
            self.assignments[target] = DiffEventAssignment()

    def set_event(self, event_hash, event_name, model_set):
        self.event = {"event_hash": event_hash, "event_name": event_name, "model_set": model_set}

    def add_trigger_species(self, species, event_hash, model_num):
        self.trigger_arrows.add({"species": species, "event_hash": event_hash}, model_num)

    def add_set_species(self, species_id, math_expr, model_num):
        self.check_target_exists(species_id)
        self.assignments[species_id].math_expr.add({"math_expr": math_expr}, model_num)

    def add_event_affect_value_arrow(self, variable_set, species, event_hash, arrow_direction, model_num):
        self.check_target_exists(variable_set)
        self.assignments[variable_set].affect_value_arrows.add(
                {"species": species, "event_hash": event_hash, "arrow_direction": arrow_direction}, model_num)

    def add_assignment_param_arrow(self, variable_set, species, event_hash, arrow_direction, model_num):
        self.check_target_exists(variable_set)
        self.assignments[variable_set].affect_value_param_arrows.add(
                {"param": species, "event_hash": event_hash, "arrow_direction": arrow_direction}, model_num)

    def add_trigger(self, math_expr, model_num):
        self.trigger_math.add({"math_expr": math_expr}, model_num)

    def add_param(self, param, event_hash, model_num):
        self.trigger_params.add({"param": param, "event_hash": event_hash}, model_num)


class DiffRule:
    def __init__(self, rule_id):
        self.rule_id = rule_id
        self.algebraic_arrows = DiffElement()
        self.modifier_arrows = DiffElement()
        self.target_arrows = DiffElement()
        self.parameter_arrows = DiffElement()
        self.rate_laws = DiffElement()

    def add_rate_law(self, model_num, converted_rate_law):
        self.rate_laws.add({"converted_rate_law": converted_rate_law}, model_num)

    def add_algebraic_arrow(self, model_num, rule_id, species_id):
        self.algebraic_arrows.add({"rule_id": rule_id, "species_id": species_id}, model_num)

    def add_modifier_arrow(self, model_num, rule_id, modifier, arrow_direction):
        self.modifier_arrows.add({"rule_id": rule_id, "modifier": modifier, "arrow_direction": arrow_direction}, model_num)

    def add_target_arrow(self, model_num, target):
        self.target_arrows.add({"target": target}, model_num)

    def add_parameter_rule(self, model_num, rule_id, param, arrow_direction):
        self.parameter_arrows.add({"rule_id": rule_id, "param": param, "arrow_direction": arrow_direction}, model_num)


class DiffReaction:
    def __init__(self):
        self.reaction = {}
        self.reactant_arrows = []
        self.product_arrows = []
        self.transcription_reaction_nodes = []
        self.transcription_product_arrows = []
        self.parameter_arrows = []

    def set_reaction(self, model_set, reaction_id, rate_law, reaction_name, converted_law,
                     fast_model_set, irreversible_model_set, product_status, is_transcription):
        self.reaction = {"model_set": model_set, "reaction_id": reaction_id, "rate_law": rate_law,
                         "reaction_name": reaction_name, "converted_law": converted_law,
                         "fast_model_set": fast_model_set, "irreversible_model_set": irreversible_model_set,
                         "product_status": product_status, "is_transcription": is_transcription}

    def add_reactant_arrow(self, model_set, reaction_id, reactant, stoich):
        self.reactant_arrows.append(
                {"model_set": model_set, "reaction_id": reaction_id, "reactant": reactant, "stoich": stoich})

    def add_product_arrow(self, model_set, reaction_id, product, stoich):
        self.product_arrows.append(
                {"model_set": model_set, "reaction_id": reaction_id, "product": product, "stoich": stoich})

    def add_transcription_reaction_node(self, model_set, reaction_id, rate_law, reaction_name, converted_law,
                                        product_status):
        self.transcription_reaction_nodes.append(
                {"model_set": model_set, "reaction_id": reaction_id, "rate_law": rate_law,
                 "reaction_name": reaction_name, "converted_law": converted_law, "product_status": product_status})

    def add_transcription_product_arrow(self, model_set, reaction_id, product, stoich):
        self.transcription_product_arrows.append(
                {"model_set": model_set, "reaction_id": reaction_id, "product": product, "stoich": stoich})

    def add_parameter_arrow(self, model_set, reaction_id, param, arrow_direction):
        self.parameter_arrows.append({"model_set": model_set, "reaction_id": reaction_id, "param": param, "arrow_direction": arrow_direction})


class DiffElement:
    def __init__(self):
        self.record = {}

    def add(self, data_tuple, model_num):
        data_tuple = FrozenDict(data_tuple)
        if data_tuple not in self.record.keys():
            self.record[data_tuple] = set()
        self.record[data_tuple].add(model_num)

    def get_models(self):
        return list(reduce(lambda x,y: x.union(self.record[y]), self.record.keys(), set()))

    def get_data(self):
        return self.record.keys()

    def all_equal(self):
        return len(self.record.values()) == 1

    def compare(self):
        if self.all_equal():
            return list(self.record.values())[0]
        else:
            return "different"

    def compare_attribute(self, attribute_name):
        val_set = False
        val = False
        for data_tuple in self.record:
            if not val_set:
                val = data_tuple[attribute_name]
                val_set = True
            elif val != data_tuple[attribute_name]:
                return "different"

        return val


class FrozenDict(collections.Mapping):
    """This class is from https://stackoverflow.com/posts/2705638/revisions"""

    def __init__(self, *args, **kwargs):
        self._d = dict(*args, **kwargs)

    def __iter__(self):
        return iter(self._d)

    def __len__(self):
        return len(self._d)

    def __getitem__(self, key):
        return self._d[key]

    def __hash__(self):
        return hash(tuple(sorted(self._d.iteritems())))