class GenerateDot():
    # This function has no dependency on BS4
    # It deals only with strings

    def __init__(self, colors, reaction_label=""):
        self.colors = colors
        self.reaction_label = reaction_label

    def assign_color(self, num_models, model_set):
        if num_models == 1:
            return "black"
        # one
        elif len(model_set) == 1 and num_models > 1:
            model_index = list(model_set)[0]
            return self.colors[model_index]
        # all
        elif len(model_set) == num_models:
            return "grey"
        # some
        elif 0 < len(model_set) < num_models:
            return "pink"

    # Used by diff_reaction()
    def print_reactant_arrow(self, num_models, model_set, reaction_id, reactant):
        color = self.assign_color(num_models, model_set)
        print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)

    def print_product_arrow(self, num_models, model_set, reaction_id, product):
        color = self.assign_color(num_models, model_set)
        print '%s -> %s [color="%s"];' % (reaction_id, product, color)

    def print_reaction_node(self, num_models, model_set, reaction_id, rate_law, reaction_name, converted_law):
        fill = ''
        if rate_law == "different":
            fill = 'fillcolor="grey", style="filled",'

        color = self.assign_color(num_models, model_set)

        if self.reaction_label == "none":
            reaction_name = ""
        elif self.reaction_label == "name":
            reaction_name = reaction_name
        elif self.reaction_label == "name+rate":
            reaction_name = reaction_name + "\n" + converted_law
        elif self.reaction_label == "rate":
            reaction_name = converted_law

        return '%s [shape="square", color="%s", %s label="%s"];' % (reaction_id, color, fill, reaction_name)

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

    def print_species_node(self, num_models, model_set, species_id, species_name):
        color = self.assign_color(num_models, model_set)
        print '"%s" [color="%s",label="%s"];' % (species_id, color, species_name)

    def print_regulatory_arrow(self, num_models, model_set, arrow_main, arrow_direction):
        color = self.assign_color(num_models, model_set)

        if arrow_direction == "monotonic_increasing":
            arrowhead = "vee"
        elif arrow_direction == "monotonic_decreasing":
            arrowhead = "tee"
        else:
            arrowhead = "dot"
        print '%s [style="dashed", color="%s", arrowhead="%s"];' % (arrow_main, color, arrowhead)
