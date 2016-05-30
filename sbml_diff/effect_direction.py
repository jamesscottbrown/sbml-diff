from bs4 import NavigableString


def collate_interactions(child_classifications):
    child_classifications = set(child_classifications)

    # all terms of same kind
    if len(child_classifications) == 1:
        return list(child_classifications)[0]

    # both kinds of monotonic terms
    if "monotonic_increasing" in child_classifications and "monotonic_decreasing" in child_classifications:
        return "mixed"

    # constant and one kind of monotonic term
    if len(child_classifications) == 2 and "monotonic_increasing" in child_classifications:
        return "monotonic_increasing"

    if len(child_classifications) == 2 and "monotonic_decreasing" in child_classifications:
        return "monotonic_decreasing"


def invert_classification(classification):

    if classification in ["mixed", "constant"]:
        return classification

    if classification == "monotonic_increasing":
        return "monotonic_decreasing"

    if classification == "monotonic_decreasing":
        return "monotonic_increasing"


def classify_basic_interaction(operator, child_classifications):
    # if any children are mixed, so is result
    if "mixed" in child_classifications:
        return "mixed"

    # monotonic function of one input:
    if operator in ["root", "exp", "ln", "log", "floor", "ceiling", "factorial", "delay"]:
        return child_classifications[0]

    if operator == "power":
        return collate_interactions(child_classifications)
        # TODO: check sign of second argument, which is in general hard

    if operator in ["plus", "times"]:
        # plus is an N-ary function
        return collate_interactions(child_classifications)

    # minus is unary or binary
    if operator == "minus":
        # unary case
        if len(child_classifications) == 1:
            invert_classification(child_classifications[0])

        # binary case
        child_classifications[1] = invert_classification(child_classifications[1])
        return collate_interactions(child_classifications)

    if operator == "divide":
        # binary operator
        child_classifications[1] = invert_classification(child_classifications[1])
        return collate_interactions(child_classifications)


def categorise_interaction(kinetic_law, species_id):
    for math in kinetic_law.select_one("math"):
        if isinstance(math, NavigableString):
            continue
        return categorise_interaction_inner(math, species_id)


def categorise_interaction_inner(expression, species_id):
    # We implicitly assume constants and powers are positive.

    # Stuff we still need to handle:
    # pi, infinity, exponential2
    # piecewise functions

    if expression.name == "cn":
        return "constant"

    if expression.name == "ci":
        if expression.string.strip() == species_id:
            return "monotonic_increasing"
        else:
            return "constant"

    # from a BS4 object, find the identity of the operator and array of it children
    operator = None
    args = []
    for child in expression.children:
        if not operator:
            operator = child.name
            if child.name == "csymbol" and child.string.strip() == "delay":
                operator = "delay"
        else:
            if isinstance(child, NavigableString):
                continue
            args.append(child)

    # classify each of the children
    child_classifications = []
    for arg in args:
        child_classifications.append(categorise_interaction_inner(arg, species_id))

    # combine these classifications based on identity of function
    return classify_basic_interaction(operator, child_classifications)
