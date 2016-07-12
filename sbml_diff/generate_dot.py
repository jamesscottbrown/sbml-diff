class GenerateDot():
    # This function has no dependency on BS4
    # It deals only with strings

    def __init__(self, colors, num_models, reaction_label="", selected_model="", show_stoichiometry=False):
        self.colors = colors
        self.num_models = num_models

        # Categorical 12-step scheme from
        # http://geog.uoregon.edu/datagraphics/color_scales.htm#Categorical%20Color%20Schemes
        #self.colors = ["#FFBF7F", "#FF7F00", "#FFFF99", "#FFFF32", "#B2FF8C", "#32FF00",
        #               "#A5EDFF", "#19B2FF", "#CCBFFF", "#654CFF", "#FF99BF", "#E51932"]

        self.selected_model = ""
        if selected_model != "":
            self.selected_model = int(selected_model) - 1

        self.show_stoichiometry = show_stoichiometry
        self.reaction_label = reaction_label

    def assign_color(self, model_set):
        if self.num_models == 1:
            return "black"
        # one
        elif len(model_set) == 1 and self.num_models > 1:
            model_index = list(model_set)[0]
            return self.colors[model_index]
        # all
        elif len(model_set) == self.num_models:
            return "grey"
        # some
        elif 0 < len(model_set) < self.num_models:
            return "black"

    def check_style(self, model_set):
        style = ', style=""'
        if self.selected_model == "" or self.selected_model in model_set:
            if len(model_set) < self.num_models:
                style = ', style="bold"'
        else:
            style = ', style="invis"'
        return style

    # Used by diff_reaction()
    def print_reactant_arrow(self, model_set, reaction_id, reactant, stoich):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        stoich_string = ''
        if self.show_stoichiometry:
            stoich_string = ', headlabel="%s", labelfontcolor=red' % stoich

        print '%s -> %s [color="%s"%s%s];' % (reactant, reaction_id, color, stoich_string, style)

    def print_product_arrow(self, model_set, reaction_id, product, stoich):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        stoich_string = ''
        if self.show_stoichiometry:
            stoich_string = ', taillabel="%s", labelfontcolor=red' % stoich

        print '%s -> %s [color="%s"%s%s];' % (reaction_id, product, color, stoich_string, style)

    def print_reaction_node(self, model_set, reaction_id, rate_law, reaction_name, converted_law):
        fill = ''
        if rate_law == "different":
            fill = 'fillcolor="grey", style="filled",'

        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        if self.reaction_label == "none":
            reaction_name = ""
        elif self.reaction_label == "name" or self.reaction_label == "":
            reaction_name = reaction_name
        elif self.reaction_label == "name+rate":
            reaction_name = reaction_name + "\n" + converted_law
        elif self.reaction_label == "rate":
            reaction_name = converted_law

        return '%s [shape="square", color="%s", %s label="%s" %s];' % (reaction_id, color, fill, reaction_name, style)

    # Used by diff_models()
    def print_header(self):
        print "\n\n"
        print "digraph comparison {"

    def print_footer(self):
        print "}"

    # Used by diff_compartment():
    def print_compartment_header(self, compartment_id):
        print "\n"
        print "subgraph cluster_%s {" % compartment_id
        print "graph[style=dotted];"
        print 'label="%s";' % compartment_id

    def print_compartment_footer(self):
        print "\n"
        print "}"

    def print_species_node(self, model_set, species_id, species_name):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print '"%s" [color="%s",label="%s" %s];' % (species_id, color, species_name, style)

    def print_regulatory_arrow(self, model_set, arrow_main, arrow_direction):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        style = style[:-1] + ',dashed"'

        if arrow_direction == "monotonic_increasing":
            arrowhead = "vee"
        elif arrow_direction == "monotonic_decreasing":
            arrowhead = "tee"
        else:
            arrowhead = "dot"
        print '%s [color="%s", arrowhead="%s" %s];' % (arrow_main, color, arrowhead, style)

    #
    def print_rule_modifier_arrow(self, model_set, rule_id, modifier):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print '%s -> rule_%s [color="%s", style="dotted" %s];' % (modifier, rule_id, color, style)

    def print_rule_target_arrow(self, model_set, rule_id, target):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print 'rule_%s -> %s [color="%s", style="dotted" %s];' % (rule_id, target, color, style)

    def print_rule_node(self, model_set, rule_id, rate_law, reaction_name, converted_rate_law):
        fill = ''
        if rate_law == "different":
            fill = 'fillcolor="grey", style="filled",'

        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        if self.reaction_label in ["none", "name", ""]:
            rule_name = ""
        elif self.reaction_label in ["name+rate", "rate"]:
            rule_name = converted_rate_law

        return 'rule_%s [shape="parallelogram", color="%s", %s label="%s" %s];' % (rule_id, color, fill, rule_name, style)

    def print_abstracted_arrow(self, model_set, modifier, target, effect_type):

        if len(model_set) == 0:
            return

        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        if effect_type in ["decrease-degredation", "increase-production"]:
            arrowhead = "vee"
        else:
            arrowhead = "tee"

        if effect_type in ["increase-degredation", "decrease-degredation"]:
            style = style[:-1] + 'dashed"'

        print '%s -> %s [style="dashed", color="%s", arrowhead="%s" %s];' % (modifier, target, color, arrowhead, style)
