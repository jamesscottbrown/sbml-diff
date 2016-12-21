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

    def generate_dot(self, diff_object):
        self.print_header()

        params_to_draw = []
        for p in diff_object.param_nodes:
            params_to_draw.append(p["variable_id"])

        for compartment_id in diff_object.compartments.keys():

            compartment = diff_object.compartments[compartment_id]

            if compartment_id is not "NONE":
                self.print_compartment_header(compartment_id)

            for species in compartment["species"]:
                self.print_species_node(species["model_set"], species["isBoundary"], species["species_id"], species["species_name"])

            for reaction in compartment["reactions"]:
                # reaction node
                r = reaction.reaction
                if r["is_transcription"]:
                    self.print_transcription_reaction_node(r["model_set"], r["reaction_id"], r["rate_law"], r["reaction_name"], r["converted_law"], r["product_status"])
                else:
                    self.print_reaction_node(r["model_set"], r["reaction_id"], r["rate_law"], r["reaction_name"], r["converted_law"],
                            r["fast_model_set"], r["irreversible_model_set"])

                # reactant arrows
                for r in reaction.reactant_arrows:
                    self.print_reactant_arrow(r["model_set"], r["reaction_id"], r["reactant"], r["stoich"])

                # parameter arrows
                for r in reaction.parameter_arrows:
                    if r["param"] in params_to_draw:
                        self.print_reaction_parameter_arrow(r["model_set"], r["reaction_id"], r["param"])

                # product arrows
                for product_arrow in reaction.product_arrows:
                    self.print_product_arrow(product_arrow["model_set"], product_arrow["reaction_id"], product_arrow["product"], product_arrow["stoich"])

                for product_arrow in reaction.transcription_product_arrows:
                    self.print_transcription_product_arrow(product_arrow["model_set"], product_arrow["reaction_id"], product_arrow["product"], product_arrow["stoich"])

            for r in compartment["regulatory_arrows"]:
                self.print_regulatory_arrow(r["model_set"], r["arrow_source"], r["arrow_target"], r["arrow_direction"])

            for r in compartment["rules"]:

                rule = r.rule
                self.print_rule_node(rule["model_set"], rule["rule_id"], rule["rate_law"], rule["converted_rate_law"])

                for arrow in r.algebraic_arrows:
                    self.print_algebraic_rule_arrow(arrow["model_set"], arrow["rule_id"], arrow["species_id"])

                for arrow in r.modifier_arrows:
                    self.print_rule_modifier_arrow(arrow["model_set"], arrow["rule_id"], arrow["modifier"], arrow["arrow_direction"])

                for arrow in r.target_arrows:
                    self.print_rule_target_arrow(arrow["model_set"], arrow["target"])

                for arrow in r.parameter_arrows:
                    if arrow["param"] in params_to_draw:
                        self.print_rule_parameter_arrow(arrow["model_set"], arrow["rule_id"], arrow["param"], arrow["arrow_direction"])

            if compartment_id is not "NONE":
                self.print_compartment_footer()

        for event in diff_object.events:
            self.print_event_diff(event, params_to_draw)

        # modified params
        for param_node in diff_object.param_nodes:
            self.print_param_node(param_node["variable_id"], param_node["variable_name"], param_node["model_set"])

        self.print_footer()

    def print_event_diff(self, event, params_to_draw):

        for r in event.trigger_params:
            if r["param"] in params_to_draw:
                self.print_event_trigger_species_arrows(r["param"], r["event_hash"], r["model_set"])

        if len(event.assignments) < 2:
            self.print_event_node(event.event["event_hash"], event.event["event_name"], event.trigger_math, event.event["model_set"])

            for s in event.trigger_arrows:
                self.print_event_trigger_species_arrows(s["species"], s["event_hash"], s["model_set"])

            for target_id in event.assignments:
                self.print_event_set_species_arrow(target_id, event.event["event_hash"], event.assignments[target_id]["model_set"])

                for a2 in event.assignments[target_id].affect_value_arrows:
                    self.print_event_affect_value_arrow(a2["species"], a2["event_hash"], a2["arrow_direction"], a2["model_set"])

        else:

            print "subgraph cluster_event_%s {\n" % event.event["event_hash"]

            # trigger arrows go to diamond event node, as in single-assignment case
            self.print_event_node(event.event["event_hash"], event.event["event_name"], event.trigger_math, event.event["model_set"])

            for s in event.trigger_arrows:
                self.print_event_trigger_species_arrows(s["species"], s["event_hash"], s["model_set"])

            for target_id in event.assignments.keys():
                # draw assignment node, like a rule
                s = event.assignments[target_id]
                rule_id = event.event["event_hash"] + "_" + target_id

                self.print_rule_node(s.model_set, rule_id, s.math_expr, s.math_expr)
                self.print_event_target_arrow(s.model_set, rule_id, target_id)

                for modifier in s.affect_value_arrows:
                    self.print_rule_modifier_arrow(modifier["model_set"], rule_id, modifier["species"], modifier["arrow_direction"])

                for param in s.affect_value_param_arrows:
                    if param["param"] in params_to_draw:
                        self.print_rule_parameter_arrow(param["model_set"], rule_id, param["param"], param["arrow_direction"])

            print "}"

    def assign_arrowhead(self, effect_direction):
        if effect_direction == "monotonic_increasing":
            arrowhead = "vee"
        elif effect_direction == "monotonic_decreasing":
            arrowhead = "tee"
        else:
            arrowhead = "none"
        return arrowhead

    def assign_color(self, model_set, ignore_difference=False):
        """
        Given a list of models containing some feature, determine what color that feature should be drawn (black if only
        one model is being considered, otherwise grey if present in all models, colored if in a single model, and black
        if in multiple models).

        Parameters
        ----------
        model_set : list of model numbers containing the feature
        ignore_difference : indicates that this function call does not imply the existence of differences between models

        Returns
        -------
        string specifying color

        """
        if self.num_models != len(model_set) and not ignore_difference:
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

    def print_reaction_parameter_arrow(self, model_set, reaction_id, param):
        style = self.check_style(model_set, 'dashed')
        color = self.assign_color(model_set)
        print '%s -> %s [color="%s" %s];' % (param, reaction_id, color, style)

    def print_transcription_reaction_node(self, model_set, reaction_id, rate_law, reaction_name, converted_law, product_status):
        base_style = ''
        if rate_law == "different":
            self.differences_found = True
            base_style = 'dashed'

        style = self.check_style(model_set, base_style)

        if self.reaction_label == "none":
            reaction_name = ""
        elif self.reaction_label == "name" or self.reaction_label == "":
            reaction_name = reaction_name
        elif self.reaction_label == "name+rate":
            reaction_name = reaction_name + "<br/>" + converted_law
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
        print result

    def print_transcription_product_arrow(self, model_set, reaction_id, product, stoich):
        """
        Draw arrow from reaction to product.

        Parameters
        ----------
        model_set : list of model numbers containing the feature

        reaction_id : id of the reaction

        product : id of the product

        stoich : stoichiometry of this product for this reaction

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
                            fast_model_set, irreversible_model_set):
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

        irreversible_model_set : list of indexes for models in which this reaction is irreversible


        """
        fill = ''
        base_style = ''
        if rate_law == "different":
            self.differences_found = True
            base_style = 'dashed'

        color = self.assign_color(model_set)
        style = self.check_style(model_set, base_style)

        if self.reaction_label == "none":
            reaction_name = ""
        elif self.reaction_label == "name" or self.reaction_label == "":
            reaction_name = reaction_name
        elif self.reaction_label == "name+rate":
            reaction_name = reaction_name + "<br/>" + converted_law
        elif self.reaction_label == "rate":
            reaction_name = converted_law

        reaction_name = self.reaction_details(reaction_name, irreversible_model_set, fast_model_set)

        print '%s [shape="rectangle", color="%s", %s label=%s %s];' % (reaction_id, color, fill, reaction_name, style)

    # Used by diff_models()
    def print_header(self):
        """ Print header needed for valid DOT file"""
        print "\n\n"
        print "digraph comparison {"
        print "rankdir = %s;" % self.rankdir

    def print_footer(self):
        """ Print footer needed for valid DOT file  """
        file_strings = []
        for i in range(0, len(self.model_names)):
            file_strings.append("<font color='%s'>%s</font>" % (self.assign_color([i], ignore_difference=True), self.model_names[i]))

        print 'label=<Files: %s>;' % ', '.join(file_strings)
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
            base_style = 'dashed'
        elif isBoundary.lower() == 'true':
            doubled = 'peripheries=2'

        style = self.check_style(model_set, base_style)
        print '"%s" [color="%s",label="%s" %s %s %s];' % (species_id, color, species_name, doubled, fill, style)


    def print_regulatory_arrow(self, model_set, arrow_source, arrow_target, arrow_direction):
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
        arrowhead = self.assign_arrowhead(arrow_direction)
        print '"%s" -> "%s" [color="%s", arrowhead="%s" %s];' % (arrow_source, arrow_target, color, arrowhead, style)

    def print_rule_modifier_arrow(self, model_set, rule_id, modifier, arrow_direction):
        """
        Draw arrow from species to rule

        Parameters
        ----------
        model_set : list of model numbers containing the feature
            
        rule_id : id of the rule
            
        modifier : id of the species affecting the rule

        arrow_direction : string representing kind of interaction - 'monotonic_increasing' (activation) or
        'monotonic_decreasing' (repression)
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set, 'dashed')

        arrowhead = self.assign_arrowhead(arrow_direction)
        print '%s -> rule_%s [color="%s", arrowhead="%s" %s];' % (modifier, rule_id, color, arrowhead, style)

    def print_rule_parameter_arrow(self, model_set, target, param, arrow_direction):
        color = self.assign_color(model_set)
        style = self.check_style(model_set, 'dashed')
        arrowhead = self.assign_arrowhead(arrow_direction)

        print '%s -> rule_%s [color="%s", arrowhead="%s" %s];' % (param, target, color, arrowhead, style)

    def print_event_target_arrow(self, model_set, event_id, target):
        """
        Draw arrow from rule to species

        Parameters
        ----------
        model_set : list of model numbers containing the feature

        target : id of the species affected by the rule
        """
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print 'rule_%s -> "%s" [color="%s", style="dotted" %s];' % (event_id, target, color, style)


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

        print 'rule_%s [shape="parallelogram", color="%s", %s label="%s" %s];' % (rule_id, color, fill, rule_name, style)

    def print_algebraic_rule_arrow(self, model_set, rule_id, species_id):
        color = self.assign_color(model_set)
        style = self.check_style(model_set)
        print 'rule_%s -> %s [color="%s", dir="none" %s];' % (rule_id, species_id, color, style)


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

    def print_event_node(self, event_hash, event_name, rate_law,  model_set):

        base_style = ''
        if rate_law == "different":
            self.differences_found = True
            base_style = 'dashed'

        if self.reaction_label == "none":
            event_name = ""
        elif self.reaction_label == "name+rate":
            event_name = event_name + "\n" + rate_law
        elif self.reaction_label == "rate":
            event_name = rate_law

        color = self.assign_color(model_set)
        style = self.check_style(model_set, base_style)

        print '"%s" [label="%s", shape="diamond", color="%s" %s];' % (event_hash, event_name, color, style)

    def print_event_trigger_species_arrows(self, species, event_hash, model_set):
        color = self.assign_color(model_set)
        print '"%s" -> "%s" [arrowhead="odot", color="%s", style="dashed"];' % (species, event_hash, color)

    def print_event_set_species_arrow(self, species_id, event_hash, model_set):
        color = self.assign_color(model_set)
        print '%s -> %s [color="%s"];' % (event_hash, species_id, color)

    def print_event_affect_value_arrow(self, species, event_hash, arrow_direction, model_set):
        color = self.assign_color(model_set)
        arrowhead = self.assign_arrowhead(arrow_direction)
        print '%s -> %s [color="%s", arrowhead="%s", style="dashed"];' % (species, event_hash, color, arrowhead)

    def print_param_node(self, variable_id, variable_name, model_set):
        color = self.assign_color(model_set)
        print '%s [label="%s", shape=none, color="%s"];' % (variable_id, variable_name, color)

    def reaction_details(self, old_label, irreversible_model_set, fast_model_set):
        """
        Add 'IR' and 'F' to label of reaction node to indicate the reaction is irreversible or fast, respectively.
        This markers are coloured independently of the rest of the node, following the same rules as other elements.

        Parameters
        ----------
        old_label : the label for the reaction (name or id, perhaps with rate expression)
        irreversible_model_set : list of indexes for models in which this reaction is irreversible
        fast_model_set : list of indexes for models in which this reaction is fast

        Returns
        -------

        """

        reversible_string = ''
        if len(irreversible_model_set) > 0:
            reversible_color = self.assign_color(irreversible_model_set)
            reversible_string = "<font color='%s'>IR</font>" % reversible_color

        fast_string = ''
        if len(fast_model_set) > 0:
            fast_color = self.assign_color(fast_model_set)
            fast_string = "<font color='%s'>F</font>" % fast_color

        if fast_string and reversible_string:
            format_string = '<%s<br/>(%s,%s)>' % (old_label, reversible_string, fast_string)
        elif fast_string:
            format_string = '<%s<br/>(%s)>' % (old_label, fast_string)
        elif reversible_string:
            format_string = '<%s<br/>(%s)>' % (old_label, reversible_string)
        else:
            format_string = '<%s>' % old_label

        if not format_string:
            format_string = '""'

        return format_string