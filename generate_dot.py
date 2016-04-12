# This function has no dependency on BS4
# It deals only with strings


def assign_color(num_models, model_set, colors):
    if num_models == 1:
        return "black"
    # one
    elif len(model_set) == 1 and num_models > 1:
        model_index = list(model_set)[0]
        return colors[model_index]
    # all
    elif len(model_set) == num_models:
        return "grey"
    # some
    elif 0 < len(model_set) < num_models:
        return "pink"


# Used by diff_reaction()
def print_reactant_arrow(num_models, model_set, colors, reactant, reaction_id):
    color = assign_color(num_models, model_set, colors)
    print '%s -> %s [color="%s"];' % (reactant, reaction_id, color)


def print_product_arrow(num_models, model_set, colors, reaction_id, product):
    color = assign_color(num_models, model_set, colors)
    print '%s -> %s [color="%s"];' % (reaction_id, product, color)


def print_rate_law(num_models, model_set, colors, rate_law, reaction_id, reaction_name):
    fill = ''
    if rate_law == "different":
        fill = 'fillcolor="grey", style="filled",'

    color = assign_color(num_models, model_set, colors)
    return '%s [shape="square", color="%s", %s label="%s"];' % (reaction_id, color, fill, reaction_name)


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
