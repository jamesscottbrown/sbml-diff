from bs4 import BeautifulSoup, NavigableString


def convert_rate_law(math, variables_not_to_substitute=False):
    return convert_rate_law_inner(math, variables_not_to_substitute)[1]


def add_parens(term_elementary, terms):
    if not term_elementary[0]:
        terms[0] = "(%s)" % terms[0]
    if not term_elementary[1]:
        terms[1] = "(%s)" % terms[1]
    return terms


def convert_rate_law_inner(expression, variables_not_to_substitute=False):

    # Stuff we still need to handle:
    # pi, infinity, exponential2
    # piecewise functions

    elementary = False
    if expression.name in ["cn", "ci"]:
        elementary = True
        term = expression.string.strip()

        if variables_not_to_substitute and term not in variables_not_to_substitute:
            term = '1.0'

        return elementary, term

    # math may contain either an <apply> or a <cn>
    if expression.name == "math":
        for child in expression.children:
            if not isinstance(child, NavigableString):
                return convert_rate_law_inner(child, variables_not_to_substitute)

    # First child is operator; next are arguments
    if expression.name == "apply":
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

        children_converted = []
        children_elementary = []
        for arg in args:
            child_elementary, child_converted = convert_rate_law_inner(arg, variables_not_to_substitute)
            children_converted.append(child_converted)
            children_elementary.append(child_elementary)

        if operator == "plus":
            return elementary, " + ".join(children_converted)
        elif operator == "minus":
            return elementary, " - ".join(children_converted)
        elif operator == "times":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, " * ".join(children_converted)
        elif operator == "divide":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s / %s" % (children_converted[0], children_converted[1])
        elif operator == "power":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s ^ %s " % (children_converted[0], children_converted[1])
        elif operator == "delay":
            return elementary, "delay(%s, %s)" % (children_converted[0], children_converted[1])
        elif operator in ["root", "exp", "ln", "log", "floor", "ceiling", "factorial"]:
            return elementary, "%s(%s)" % (operator, children_converted[0])
