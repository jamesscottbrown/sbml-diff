from bs4 import BeautifulSoup, NavigableString

def convert_rate_law(math):
    return convert_rate_law_inner(math)[1]


def add_parens(term_elementary, terms):
    if not term_elementary[0]:
        terms[0] = "(%s)" % terms[0]
    if not term_elementary[1]:
        terms[1] = "(%s)" % terms[1]
    return terms


def convert_rate_law_inner(expression):

    # Stuff we still need to handle:
    # pi, infinity, exponential2
    # delay csymbol
    # piecewise functions

    elementary = False
    if expression.name in ["cn", "ci"]:
        elementary = True
        return elementary, expression.string.strip()

    # math may contain either an <apply> or a <cn>
    if expression.name == "math":
        for child in expression.children:
            if not isinstance(child, NavigableString):
                return convert_rate_law_inner(child)

    # First child is operator; next are arguments
    if expression.name == "apply":
        operator = None
        args = []
        for child in expression.children:
            if not operator:
                operator = child.name
            else:
                if isinstance(child, NavigableString):
                    continue
                args.append(child)

        children_converted = []
        children_elementary = []
        for arg in args:
            child_elementary, child_converted = convert_rate_law_inner(arg)
            children_converted.append(child_converted)
            children_elementary.append(child_elementary)

        if operator == "plus":
            return elementary, " + ".join(children_converted)
        elif operator == "minus":
            return elementary, " - ".join(children_converted)
        elif operator == "times":
            return elementary, " * ".join(children_converted)
        elif operator == "divide":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s / %s" % (children_converted[0], children_converted[1])
        elif operator == "power":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s ^ %s " % (children_converted[0], children_converted[1])
        elif operator in ["root", "exp", "ln", "log", "floor", "ceiling", "factorial"]:
            return elementary, "%s(%s)" % (operator, children_converted[0])
