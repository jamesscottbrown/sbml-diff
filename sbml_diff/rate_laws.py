from bs4 import BeautifulSoup, NavigableString


def convert_rate_law(math, variables_not_to_substitute=False, executable=False):
    """
    A wrapper for convert_rate_law_inner that returns only the converted expression.

    Parameters
    ----------
    math : BeautifulSoup object representing a rateLaw (a bs4.element.Tag)

    variables_not_to_substitute : if specified, the name of any species whose id is not in this list is replaced by 1.0
         (Default value = False)

    executable : if True, return a less human-readable string that can be eval'ed in Python (e.g. containing Math.e
        instead of e) (Default value = False)

    Returns
    -------
    string representation of the kineticLaw

    """
    return convert_rate_law_inner(math, variables_not_to_substitute, executable)[1]


def add_parens(term_elementary, terms):
    """
    If any elements in the first argument is false, wrap the corresponding elements of the second argument in parentheses.

    Parameters
    ----------
    term_elementary : boolean - if this is false, wrap terms in parentheses
        
    terms : string

    """
    for ind, is_elementary in enumerate(term_elementary):
        if not is_elementary:
            terms[ind] = "(%s)" % terms[ind]
    return terms


def convert_rate_law_inner(expression, variables_not_to_substitute=False, executable=False):
    """
    Recursively convert a MathML expression to a string.
    Limitations: we do not handle piecewise functions or user-defined functions.

    Parameters
    ----------
    expression :
        
    variables_not_to_substitute : if specified, the name of any species whose id is not in this list is replaced by 1.0
         (Default value = False)

    executable : if True, return a less human-readable string that can be eval'ed in Python (e.g. containing Math.e
        instead of e) (Default value = False)

    Returns
    -------
    string representation of the kineticLaw

    """

    executable_replacement = {'exp': 'math.exp', 'ln': 'math.log', 'log': 'math.log10', 'ceiling': 'math.ceil',
                              'floor': 'math.floor', 'factorial': 'math.factorial', 'pi': 'math.pi', 'e': 'math.e',
                              'infinity': 'float("Inf")', "sqrt": "math.sqrt", "abs": "abs"}

    elementary = False

    if expression.name == "cn":

        if "type" in expression.attrs.keys():
            children = []
            for child in expression.children:
                children.append(child)

            if expression.attrs["type"] == "e-notation":
                term = "%s * 10^(%s)" % (children[0], children[2])
                if executable:
                    term = "%s * 10**(%s)" % (children[0], children[2])
            elif expression.attrs["type"] in ["real", "integer"]:
                term = children[0]
            elif expression.attrs["type"] == "rational":
                term = "%s/%s" % (children[0], children[2])
        else:
            term = expression.string.strip()

        elementary = True

        if variables_not_to_substitute and term not in variables_not_to_substitute:
            term = '1.0'

        return elementary, term

    elif expression.name == "ci":
        elementary = True
        term = expression.string.strip()

        if variables_not_to_substitute and term not in variables_not_to_substitute:
            term = '1.0'

        return elementary, term

    if expression.name in ["pi", "infinity"]:
        if executable:
            return True, executable_replacement[expression.name]
        return True, expression.name
    if expression.name == "exponentiale":
        if executable:
            return True, executable_replacement["e"]
        return True, "e"

    # math may contain either an <apply> or a <cn>
    if expression.name == "math":
        for child in expression.children:
            if not isinstance(child, NavigableString):
                return convert_rate_law_inner(child, variables_not_to_substitute, executable)

    if expression.name == "csymbol":
        if "time" in expression.attrs['definitionURL']:
            if executable:
                return True, '1'
            return True, "t"
        if "avogadro" in expression.attrs['definitionURL']:
            if executable:
                return True, '1'
            return True, "N_A"

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
            child_elementary, child_converted = convert_rate_law_inner(arg, variables_not_to_substitute, executable)
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
            if executable:
                return elementary, "%s ** %s " % (children_converted[0], children_converted[1])
            return elementary, "%s ^ %s " % (children_converted[0], children_converted[1])
        elif operator == "delay":
            if executable:
                return elementary, children_converted[0]
            return elementary, "delay(%s, %s)" % (children_converted[0], children_converted[1])
        elif operator in ["exp", "ln", "log", "floor", "ceiling", "factorial", "abs"]:
            if executable:
                return elementary, "%s(%s)" % (executable_replacement(operator), children_converted[0])
            return elementary, "%s(%s)" % (operator, children_converted[0])
        elif operator == "root":

            # default to sqrt()
            if len(children_converted) == 1:
                if executable:
                    return elementary, "%s(%s)" % (executable_replacement("sqrt"), children_converted[0])
                return elementary, "%s(%s)" % ("sqrt", children_converted[0])

            # otherwise root(n,a) = pow(a, 1/n)
            if len(children_converted) == 2:
                if executable:
                    return elementary, "pow(%s, 1/%s)" % (children_converted[1], children_converted[0])
                return elementary, "%s(%s, %s)" % ("root", children_converted[0], child_converted[1])
