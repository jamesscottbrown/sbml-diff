class GenerateDot:
    """This class actually generates the DOT output.
    
    It has no dependency on BeautifulSoup, and works with strings, rather than BeautifulSoup objects.

    The print_ functions accept an argument model_set, which specifies which models contain the corresponding feature.
    """

    def __init__(self, colors, num_models, reaction_label="", selected_model="", show_stoichiometry=False):
        self.colors = colors
        self.num_models = num_models

        # If too few colors specified, extend using categorical 12-step scheme from
        # http://geog.uoregon.edu/datagraphics/color_scales.htm#Categorical%20Color%20Schemes
        default_colors = ["#FFBF7F", "#FF7F00", "#FFFF99", "#FFFF32", "#B2FF8C", "#32FF00",
                       "#A5EDFF", "#19B2FF", "#CCBFFF", "#654CFF", "#FF99BF", "#E51932"]
        if len(self.colors) < self.num_models:
            spare_colors = set(default_colors).difference(self.colors)
            extra_colors = self.num_models - len(self.colors)
            self.colors.extend(spare_colors[1:extra_colors])

        self.selected_model = ""
        if selected_model != "":
            self.selected_model = int(selected_model) - 1

        self.show_stoichiometry = show_stoichiometry
        self.reaction_label = reaction_label

    def assign_color(self, model_set):
        """
        Given a list of models containing some feature, determine what color that feature should be drawn (black if only
        one model is being considered, otherwise grey if present in all models, colored if in a single model, and black
        if in multiple models).

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            

        Returns
        -------
        string specifying color

        """
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

    def check_style(self, model_set, base_style=''):
        """
        Determine whether a feature should be drawn in bold (because it is not in all models), or invisible (because the
        selected model(s) do not contain it), and add these details to the 'style' string.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        base_style : other style attributes that must be applied (e.g. dashed, or a fillcolor)
             (Default value = '')

        Returns
        -------
        a string of the form ', style="something"'

        """
        style = ', style="%s"' % base_style

        base_style = "," + base_style
        if self.selected_model == "" or self.selected_model in model_set:
            if len(model_set) < self.num_models:
                style = ', style="bold%s"' % base_style
        else:
            style = ', style="invis%s"' % base_style
        return style

    def print_reactant_arrow(self, model_set, reaction_id, reactant, stoich):
        """
        Draw arrow from reactant to reaction.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        reaction_id : id of the reaction
            
        reactant : id of the reactant
            
        stoich : stoichiometry of this reactant for this reaction


        Returns
        -------
        string representing this arrow

        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        stoich_string = ''
        if self.show_stoichiometry:
            stoich_string = ', headlabel="%s", labelfontcolor=red' % stoich

        print '%s -> %s [color="%s"%s%s];' % (reactant, reaction_id, color, stoich_string, style)

    def print_product_arrow(self, model_set, reaction_id, product, stoich):
        """
        Draw arrow from reaction to product.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        reaction_id : id of the reaction
            
        product : id of the product
            
        stoich : stoichiometry of this product for this reaction


        Returns
        -------
        string representing this arrow

        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        stoich_string = ''
        if self.show_stoichiometry:
            stoich_string = ', taillabel="%s", labelfontcolor=red' % stoich

        print '%s -> %s [color="%s"%s%s];' % (reaction_id, product, color, stoich_string, style)

    def print_reaction_node(self, model_set, reaction_id, rate_law, reaction_name, converted_law):
        """
        Draw square node representing a reaction.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        reaction_id : id of the reaction
            
        rate_law : bs4.element.Tag specifying the kineticLaw
            
        reaction_name : name of the reaction
            
        converted_law : human-readable string representation of the kineticLaw


        Returns
        -------
        string representing this node

        """
        fill = ''
        base_style = ''
        if rate_law == "different":
            fill = 'fillcolor="grey",'
            base_style = 'filled'

        color = self.assign_color(model_set)
        style = self.check_style(model_set, base_style)

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
        """ Print header needed for valid DOT file"""
        print "\n\n"
        print "digraph comparison {"

    def print_footer(self):
        """ Print footer needed for valid DOT file  """
        print "}"

    def print_compartment_header(self, compartment_id):
        """
        Print DOT code to create a new subgraph representing a compartment.

        Parameters
        ----------
        compartment_id : id of a compartment
        """
        print "\n"
        print "subgraph cluster_%s {" % compartment_id
        print "graph[style=dotted];"
        print 'label="%s";' % compartment_id

    def print_compartment_footer(self):
        """ Print DOT code to end the subgraph representing a compartment """
        print "\n"
        print "}"

    def print_species_node(self, model_set, species_id, species_name):
        """
        Draw node representing a species.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        species_id : id of a species
            
        species_name : name of a species
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print '"%s" [color="%s",label="%s" %s];' % (species_id, color, species_name, style)

    def print_regulatory_arrow(self, model_set, arrow_main, arrow_direction):
        """
        Draw arrow corresponding to regulatory interaction affecting reaction

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        arrow_main : the DOT edge_stmt for the edge (eg. 'A -> B')
            
        arrow_direction : string representing kind of interaction - 'monotonic_increasing' (activation) or
        'monotonic_decreasing' (repression)
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set, 'dashed')

        if arrow_direction == "monotonic_increasing":
            arrowhead = "vee"
        elif arrow_direction == "monotonic_decreasing":
            arrowhead = "tee"
        else:
            arrowhead = "dot"
        print '%s [color="%s", arrowhead="%s" %s];' % (arrow_main, color, arrowhead, style)

    def print_rule_modifier_arrow(self, model_set, rule_id, modifier):
        """
        Draw arrow from species to rule

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        rule_id : id of the rule
            
        modifier : id of the species affecting the rule
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print '%s -> rule_%s [color="%s", style="dotted" %s];' % (modifier, rule_id, color, style)

    def print_rule_target_arrow(self, model_set, rule_id, target):
        """
        Draw arrow from rule to species

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        rule_id : id of the rule
            
        target : id of the species affected by the rule
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print 'rule_%s -> %s [color="%s", style="dotted" %s];' % (rule_id, target, color, style)

    def print_rule_node(self, model_set, rule_id, rate_law, reaction_name, converted_rate_law):
        """
        Draw node corresponding to rule.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        rule_id : id of the rule
            
        rate_law :
            
        reaction_name :
            
        converted_rate_law :


        Returns
        -------
        string representing the node

        """
        fill = ''
        base_style = ''
        if rate_law == "different":
            fill = 'fillcolor="grey",'
            base_style = 'filled'

        color = self.assign_color(model_set)
        style = self.check_style(model_set, base_style)

        rule_name = ""
        if self.reaction_label in ["name+rate", "rate"]:
            rule_name = converted_rate_law

        return 'rule_%s [shape="parallelogram", color="%s", %s label="%s" %s];' % (rule_id, color, fill, rule_name, style)

    def print_abstracted_arrow(self, model_set, modifier, target, effect_type):
        """
        Draw an arrow between two species, indicating an interaction.
        Used to produce an abstract diagram.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        modifier : id of species affecting the target
            
        target : id of species affected by this interaction
            
        effect_type : type of interaction ("increase-degredation", "decrease-degredation",
        "decrease-degredation", "increase-production")
        """

        if len(model_set) == 0:
            return

        base_style = ''
        if effect_type in ["increase-degredation", "decrease-degredation"]:
            base_style = 'dashed'

        color = self.assign_color(model_set)
        style = self.check_style(model_set, base_style)

        if effect_type in ["decrease-degredation", "increase-production"]:
            arrowhead = "vee"
        else:
            arrowhead = "tee"

        print '%s -> %s [style="dashed", color="%s", arrowhead="%s" %s];' % (modifier, target, color, arrowhead, style)
