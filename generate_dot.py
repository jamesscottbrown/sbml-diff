# This function has no dependency on BS4
# It deals only with strings


def assign_color(num_models, model_set, colors):
    # one
    if len(model_set) == 1 and num_models > 1:
        model_index = list(model_set)[0]
        return colors[model_index]
    # all
    elif len(model_set) == num_models:
        return "grey"

    # some - hardcoded pink


# Used by diff_reaction_common()
def print_reactant_arrow(num_models, model_set, colors, reactant, reaction_id):

    # one
    if len(model_set) == 1 and num_models > 1:
        color = assign_color(num_models, model_set, colors)
        print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)

    # all
    if len(model_set) == num_models:
        print '%s -> %s [color="grey"];' % (reactant, reaction_id)

    # some
    if 0 < len(model_set) < num_models:
        print '%s -> %s [color="pink"];' % (reactant, reaction_id)


def print_product_arrow(num_models, model_set, colors, reaction_id, product):
    # one
    if len(model_set) == 1 and num_models > 1:
        color = assign_color(num_models, model_set, colors)
        print '%s -> %s [color="%s"];' % (reaction_id, product, color)

    # all
    if len(model_set) == num_models:
        print '%s -> %s [color="grey"];' % (reaction_id, product)

    # some
    if 0 < len(model_set) < num_models:
        print '%s -> %s [color="pink"];' % (reaction_id, product)


def print_rate_law(rate_law, reaction_id, reaction_name):
    if rate_law == "different":
        return '%s [shape="square", fillcolor="grey", style="filled", label="%s"];' % (reaction_id, reaction_name)
    else:
        return '%s [shape="square", color="grey", label="%s"];' % (reaction_id, reaction_name)


# Used by diff_models()
def print_header():
    print "\n\n"
    print "digraph comparison {"


def print_footer():
    print "}"


# Used by diff_compartment():
def print_compartment_header(compartment_id):
    print "\n"
    print "subgraph cluster_%s {" % compartment_id
    print "graph[style=dotted];"
    print 'label="%s";' % compartment_id


def print_compartment_footer():
    print "\n"
    print "}"


def print_species(num_models, colors, species, species_status, species_name):
    color = assign_color(num_models, species_status, colors)
    print '"%s" [color="%s",label="%s"];' % (species, color, species_name)


def print_regulatory_arrow(arrow_direction, arrow_status, arrow_main, colors, num_models):
    color = assign_color(num_models, arrow_status, colors)

    if arrow_direction == "monotonic_increasing":
        arrowhead = "vee"
    elif arrow_direction == "monotonic_decreasing":
        arrowhead = "tee"
    else:
        arrowhead = "dot"
    print '%s [style="dashed", color="%s", arrowhead="%s"];' % (arrow_main, color, arrowhead)


# Used by diff_reactions()
def print_reaction(num_models, model_set, colors, reactant_list, product_list, reaction_id, reaction_name):
    reaction_string = ""

    # one
    if len(model_set) == 1 and num_models > 1:
        color = assign_color(num_models, model_set, colors)

        for reactant in reactant_list:
            print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)

        for product in product_list:
            print '%s -> %s [color="%s"];' % (reaction_id, product, color)

        reaction_string = '%s [shape="square", color="%s", label="%s"];' % (reaction_id, color, reaction_name)

    # all - handled before this function is called

    # some
    if 1 < len(model_set) < num_models:
        reaction_string = '%s [shape="square", color="pink", label="%s"];' % (reaction_id, reaction_name)

    return reaction_string
