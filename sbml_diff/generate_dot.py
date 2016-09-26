class GenerateDot:
    """This class actually generates the DOT output.
    
    It has no dependency on BeautifulSoup, and works with strings, rather than BeautifulSoup objects.

    The print_ functions accept an argument model_set, which specifies which models contain the corresponding feature.
    """

    def __init__(self, colors, num_models, reaction_label="", selected_model="", show_stoichiometry=False, rankdir="TB",
                 model_names=False):
        """

        Parameters
        ----------
        colors : list of colors corresponding to each model (see http://graphviz.org/doc/info/colors.html)
        num_models : number of models being compared
        reaction_label : option specifying how reaction nodes are labelled ("none"/"name"/"rate"/"name+rate")
        selected_model : if this is specified, any feature that is not in this model is given style 'invis'
        show_stoichiometry : if true, arrow between species and reaction nodes are labelled with stoichiometric coefficient

        """
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

        if not model_names:
            model_names = [""] * num_models
        self.model_names = model_names

        self.show_stoichiometry = show_stoichiometry
        self.reaction_label = reaction_label
        self.rankdir = rankdir
        self.differences_found = False

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
        if self.num_models != len(model_set):
            self.differences_found = True

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
        if self.show_stoichiometry or stoich == '?':
            stoich_string = ', headlabel="%s", labelfontcolor=red' % stoich

        if stoich == '?':
            color = "black"
            self.differences_found = True

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
        if self.show_stoichiometry or stoich == '?':
            stoich_string = ', taillabel="%s", labelfontcolor=red' % stoich

        if stoich == '?':
            color = "black"
            self.differences_found = True

        print '%s -> %s [color="%s"%s%s];' % (reaction_id, product, color, stoich_string, style)

    def print_transcription_reaction_node(self, model_set, reaction_id, rate_law, reaction_name, converted_law, product_status):
        base_style = ''
        if rate_law == "different":
            self.differences_found = True
            fill = 'fillcolor="grey",'
            base_style = 'filled'

        style = self.check_style(model_set, base_style)

        if self.reaction_label == "none":
            reaction_name = ""
        elif self.reaction_label == "name" or self.reaction_label == "":
            reaction_name = reaction_name
        elif self.reaction_label == "name+rate":
            reaction_name = reaction_name + "\n" + converted_law
        elif self.reaction_label == "rate":
            reaction_name = converted_law

        products = product_status.keys()
        result = ""
        result += "subgraph cluster_%s {\n" % reaction_id
        result += 'label="%s";\n' % reaction_name
        result += style[1:] + ";\n"

        result += 'color="%s";\n' % self.assign_color(model_set)

        for product in product_status:
            color = self.assign_color(product_status[product])
            result += 'cds_%s_%s [fillcolor="%s", style=filled, color="black", shape="cds", label=""];\n' % (reaction_id, product, color)

        result += '%s [shape=promoter, label=""];\n' % reaction_id
        result += '%s -> cds_%s_%s [arrowhead="none"];\n' % (reaction_id, reaction_id, products[0])
        for i in range(len(products)-1):
            result += "cds_%s_%s -> cds_%s_%s;\n" % (reaction_id, products[i], reaction_id, products[i-1])

        result += "}\n\n"
        return result

    def print_transcription_product_arrow(self, model_set, reaction_id, product, stoich):
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
        if self.show_stoichiometry or stoich == '?':
            stoich_string = ', taillabel="%s", labelfontcolor=red' % stoich

        if stoich == '?':
            color = "black"
            self.differences_found = True

        print 'cds_%s_%s -> %s [color="%s"%s%s];' % (reaction_id, product, product, color, stoich_string, style)

    def print_reaction_node(self, model_set, reaction_id, rate_law, reaction_name, converted_law,
                            fast_model_set, reversible_model_set):
        """
        Draw rectangular node representing a reaction.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        reaction_id : id of the reaction
            
        rate_law : bs4.element.Tag specifying the kineticLaw
            
        reaction_name : name of the reaction
            
        converted_law : human-readable string representation of the kineticLaw

        fast_model_set : list of indexes for models in which this reaction is fast

        reversible_model_set : list of indexes for models in which this reaction is reversible

        Returns
        -------
        string representing this node

        """
        fill = ''
        base_style = ''
        if rate_law == "different":
            self.differences_found = True
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

        reaction_name = self.reaction_details(reaction_name, reversible_model_set, fast_model_set)

        return '%s [shape="rectangle", color="%s", %s label=%s %s];' % (reaction_id, color, fill, reaction_name, style)

    # Used by diff_models()
    def print_header(self):
        """ Print header needed for valid DOT file"""
        print "\n\n"
        print "digraph comparison {"
        print "rankdir = %s;" % self.rankdir

    def print_footer(self):
        """ Print footer needed for valid DOT file  """
        print 'label="Files: %s";' % ', '.join(self.model_names)
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

    def print_species_node(self, model_set, isBoundary, species_id, species_name):
        """
        Draw node representing a species.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        species_id : id of a species
            
        species_name : name of a species
        """
        color = self.assign_color(model_set)

        fill = ""
        base_style = ""
        doubled = ""

        if isBoundary == '?':
            fill = ', fillcolor="grey"'
            base_style = 'filled'
        elif isBoundary.lower() == 'true':
            doubled = 'peripheries=2'

        style = self.check_style(model_set, base_style)
        print '"%s" [color="%s",label="%s" %s %s %s];' % (species_id, color, species_name, doubled, fill, style)


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
        elif arrow_direction == "?":
            arrowhead = "odiamond"
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

    def print_rule_target_arrow(self, model_set, target):
        """
        Draw arrow from rule to species

        Parameters
        ----------
        model_set : list of model numbers containing the feature

        target : id of the species affected by the rule
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print 'rule_%s -> %s [color="%s", style="dotted" %s];' % (target, target, color, style)

    def print_rule_node(self, model_set, rule_id, rate_law, converted_rate_law):
        """
        Draw node corresponding to rule.

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        rule_id : id of the rule
            
        rate_law :

        converted_rate_law :


        Returns
        -------
        string representing the node

        """
        fill = ''
        base_style = ''
        if rate_law == "different":
            self.differences_found = True
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

    def print_event_node(self, event_hash, event_name, model_set):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)

        if self.reaction_label == "none":
            event_name = ""

        print '%s [label="%s", shape="diamond", color="%s" %s];' % (event_hash, event_name, color, style)

    def print_event_trigger_species_arrows(self, species, event_hash, model_set):
        color = self.assign_color(model_set)
        print '%s -> %s [arrowhead="diamond", color="%s"];' % (species, event_hash, color)

    def print_event_set_species_arrow(self, species_id, event_hash, model_set):
        color = self.assign_color(model_set)
        print '%s -> %s [color="%s"];' % (event_hash, species_id, color)

    def print_event_affect_value_arrow(self, species, event_hash, model_set):
        color = self.assign_color(model_set)
        print '%s -> %s [color="%s", style="dashed"];' % (species, event_hash, color)

    def print_param_node(self, variable_id, variable_name, model_set):
        color = self.assign_color(model_set)
        print '%s [label="%s", shape=none, color=%s];' % (variable_id, variable_name, color)

    def reaction_details(self, old_label, reversible_model_set, fast_model_set):
        """
        Add 'R' and 'F' to label of reaction node to indicate the reaction is reversible or fast, respectively.
        This markers are coloured independently of the rest of the node, following the same rules as other elements.

        Parameters
        ----------
        old_label : the label for the reaction (name or id, perhaps with rate expression)
        reversible_model_set : list of indexes for models in which this reaction is reversible
        fast_model_set : list of indexes for models in which this reaction is fast

        Returns
        -------

        """

        reversible_string = ''
        if len(reversible_model_set) > 0:
            reversible_color = self.assign_color(reversible_model_set)
            reversible_string = '<<font color="%s">R</font>' % reversible_color

        fast_string = ''
        if len(fast_model_set) > 0:
            fast_color = self.assign_color(fast_model_set)
            fast_string = "<font color='%s'>F</font>" % fast_color

        if fast_string and reversible_string:
            format_string = '<%s<br/>%s,%s>' % (old_label, reversible_string, fast_string)
        elif fast_string:
            format_string = '<%s<br/>%s>' % (old_label, fast_string)
        elif reversible_string:
            format_string = '<%s<br/>%s>' % (old_label, reversible_string)
        else:
            format_string = old_label

        return format_string