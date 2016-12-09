class DiffObject:
    def __init__(self):
        self.compartments = {}
        self.add_compartment("NONE")
        self.events = []  # should probably be moved into compartment ?
        self.rules = []
        self.reactions = []
        self.param_nodes = []

    def add_compartment(self, compartment):
        self.compartments[compartment] = {}
        self.compartments[compartment]["rules"] = []
        self.compartments[compartment]["reactions"] = []
        self.compartments[compartment]["species"] = []
        self.compartments[compartment]["regulatory_arrows"] = []

    def check_compartment_exists(self, compartment):
        if compartment not in self.compartments.keys():
            self.add_compartment(compartment)

    def add_event(self):
        new_event = DiffEvent()
        self.events.append(new_event)
        return new_event

    def add_rule(self, compartment=""):
        new_rule = DiffRule()

        if not compartment:
            compartment = "NONE"

        self.check_compartment_exists(compartment)
        self.compartments[compartment]["rules"].append(new_rule)
        return new_rule

    def add_reaction(self, compartment=""):
        new_reaction = DiffReaction()

        if not compartment:
            compartment = "NONE"
        self.check_compartment_exists(compartment)
        self.compartments[compartment]["reactions"].append(new_reaction)
        return new_reaction

    def add_species(self, compartment, model_set, is_boundary, species_id, species_name):
        self.check_compartment_exists(compartment)
        self.compartments[compartment]["species"].append(
                {"model_set": model_set, "isBoundary": is_boundary, "species_id": species_id,
                 "species_name": species_name})

    def add_regulatory_arrow(self, compartment, model_set, arrow_source, arrow_target, arrow_direction):
        self.check_compartment_exists(compartment)
        self.compartments[compartment]["regulatory_arrows"].append(
                {"model_set": model_set, "arrow_source": arrow_source, "arrow_target": arrow_target,
                 "arrow_direction": arrow_direction})

    def add_param_node(self, variable_id, variable_name, model_set):
        self.param_nodes.append({"variable_id": variable_id, "variable_name": variable_name, "model_set": model_set})


class DiffEvent:
    def __init__(self):
        self.event = {}
        self.trigger_arrows = []
        self.set_species_arrows = []
        self.affect_value_arrows = []

    def set_event(self, event_hash, event_name, model_set):
        self.event = {"event_hash": event_hash, "event_name": event_name, "model_set": model_set}

    def add_trigger_species(self, species, event_hash, model_set):
        self.trigger_arrows.append({"species": species, "event_hash": event_hash, "model_set": model_set})

    def add_set_species(self, species_id, event_hash, model_set):
        self.set_species_arrows.append({"species_id": species_id, "event_hash": event_hash, "model_set": model_set})

    def add_event_affect_value_arrow(self, species, event_hash, arrow_direction, model_set):
        self.affect_value_arrows.append(
                {"species": species, "event_hash": event_hash, "arrow_direction": arrow_direction,
                 "model_set": model_set})


class DiffRule:
    def __init__(self):
        self.rule = {}
        self.assingment_arrows = []
        self.modifier_arrows = []
        self.target_arrows = []

    def set_rule(self, model_set, rule_id, rate_law, converted_rate_law):
        self.rule = {"model_set": model_set, "rule_id": rule_id, "rate_law": rate_law,
                     "converted_rate_law": converted_rate_law}

    def add_assignment_arrow(self, model_set, rule_id, species_id):
        self.assingment_arrows.append({"model_set": model_set, "rule_id": rule_id, "species_id": species_id})

    def add_modifier_arrow(self, model_set, rule_id, modifier, arrow_direction):
        self.modifier_arrows.append(
                {"model_set": model_set, "rule_id": rule_id, "modifier": modifier, "arrow_direction": arrow_direction})

    def add_target_arrow(self, model_set, target):
        self.target_arrows.append({"model_set": model_set, "target": target})


class DiffReaction:
    def __init__(self):
        self.reaction = {}
        self.reactant_arrows = []
        self.product_arrows = []
        self.transcription_reaction_nodes = []
        self.transcription_product_arrows = []

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
