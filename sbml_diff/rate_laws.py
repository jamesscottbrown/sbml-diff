from bs4 import BeautifulSoup, NavigableString
import copy
import sys


def convert_rate_law(math, initial_values=False, non_default_variables=False, non_default_values=1, output_type=""):
    """
    A wrapper for convert_rate_law_inner that returns only the converted expression.

    Parameters
    ----------
    math : BeautifulSoup object representing a rateLaw (a bs4.element.Tag)

    non_default_variables : if specified, the name of any species whose id is not in this list is replaced by 1.0
         (Default value = False)

    output_type : if "executable", return a less human-readable string that can be eval'ed in Python; if "sympy"
        generate string that uses sympy functions (rather than math functions)

    Returns
    -------
    string representation of the kineticLaw

    """
    if not initial_values:
        initial_values = {}

    if not math:
        return ""

    if math.select_one('piecewise') or math.name == 'piecewise':
        sys.stderr.write("Encountered a piecewise function\n")
        if output_type in ["executable", "sympy"]:
            return "piecewise"
        return ""

    return convert_rate_law_inner(math, initial_values, non_default_variables, non_default_values, output_type)[1]


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


def convert_function(output_type, function_name):
    executable_replacement = {'exp': 'math.exp', 'ln': 'math.log', 'log': 'math.log10', 'ceiling': 'math.ceil',
                          'floor': 'math.floor', 'factorial': 'math.factorial', 'pi': 'math.pi', 'e': 'math.e',
                          'infinity': 'float("Inf")', "sqrt": "math.sqrt", "abs": "abs", "cos": "math.cos",
                          "sin": "math.sin", "tan": "math.tan", "sinh": "math.sinh", "cosh": "math.cosh",
                          "tanh": "math.tanh", "arcsin": "math.asin", "arccos": "math.acos", "arctan": "math.atan",
                          "t": "1", "N_A": "1"}

    sympy_replacement = {'exp': 'sympy.exp', 'ln': 'sympy.log', 'log': 'sympy.log10', 'ceiling': 'sympy.ceiling',
                         'floor': 'sympy.floor', 'factorial': 'sympy.factorial', 'pi': 'sympy.pi', 'e': 'sympy.mpmath.e',
                         'infinity': 'float("Inf")', "sqrt": "sympy.sqrt", "abs": "abs", "cos": "sympy.cos",
                         "sin": "sympy.sin", "tan": "sympy.tan", "sinh": "sympy.sinh", "cosh": "sympy.cosh",
                         "tanh": "sympy.tanh", "arcsin": "sympy.asin", "arccos": "sympy.acos", "arctan": "sympy.atan",
                         "root": "sympy.root",
                         "t": "t", "N_A": "N_A"}

    if not output_type:
        return function_name
    elif output_type == "executable":
        return executable_replacement[function_name]
    elif output_type == "sympy":
        return sympy_replacement[function_name]


def convert_rate_law_inner(expression, initial_values, non_default_variables=False, non_default_values=1, output_type=""):
    """
    Recursively convert a MathML expression to a string.
    Limitations: we do not handle piecewise functions or user-defined functions.

    Parameters
    ----------
    expression :
        
    non_default_variables : if specified, the name of any species whose id is not in this list is replaced by 1.0
         (Default value = False)

    executable : if True, return a less human-readable string that can be eval'ed in Python (e.g. containing Math.e
        instead of e) (Default value = False)

    Returns
    -------
    string representation of the kineticLaw

    """

    elementary = False

    generate_code = (output_type in ["executable", "sympy"])

    if expression.name == "cn":

        if "type" in list(expression.attrs.keys()):
            children = []
            term = ""
            for child in expression.children:
                children.append(child)

            if expression.attrs["type"] == "e-notation":
                term = "%s * 10^(%s)" % (children[0], children[2])
                if generate_code:
                    term = "%s * 10**(%s)" % (children[0], children[2])
            elif expression.attrs["type"] in ["real", "integer"]:
                term = children[0]
            elif expression.attrs["type"] == "rational":
                term = "%s/%s" % (children[0], children[2])
        else:
            term = expression.string.strip()

        elementary = True

        return elementary, term

    elif expression.name == "ci":
        elementary = True
        term = expression.string.strip()

        if non_default_variables:
            if term in non_default_variables:
                term = non_default_values
            elif term in list(initial_values.keys()):
                term = initial_values[term]
            else:
                term = '1.0'

        return elementary, term

    if expression.name in ["pi", "infinity"]:
        return True, convert_function(output_type, expression.name)
    if expression.name == "exponentiale":
        return True, convert_function(output_type, "e")

    # math may contain either an <apply> or a <cn>
    if expression.name == "math":
        for child in expression.children:
            if not isinstance(child, NavigableString):
                return convert_rate_law_inner(child, initial_values, non_default_variables, non_default_values, output_type)

    if expression.name == "csymbol":
        if "time" in expression.attrs['definitionURL']:
            return True, convert_function(output_type, "t")

        if "avogadro" in expression.attrs['definitionURL']:
            return True, convert_function(output_type, "N_A")

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
            child_elementary, child_converted = convert_rate_law_inner(arg, initial_values, non_default_variables, non_default_values, output_type)
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
            if generate_code:
                return elementary, "%s ** %s " % (children_converted[0], children_converted[1])
            return elementary, "%s ^ %s " % (children_converted[0], children_converted[1])
        elif operator == "gt":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s > %s " % (children_converted[0], children_converted[1])
        elif operator == "geq":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s >= %s " % (children_converted[0], children_converted[1])
        elif operator == "eq":
            children_converted = add_parens(children_elementary, children_converted)
            if generate_code:
                return elementary, "%s == %s " % (children_converted[0], children_converted[1])
            return elementary, "%s = %s " % (children_converted[0], children_converted[1])
        elif operator == "lt":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s < %s " % (children_converted[0], children_converted[1])
        elif operator == "leq":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s <= %s " % (children_converted[0], children_converted[1])
        elif operator == "or":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s or %s " % (children_converted[0], children_converted[1])
        elif operator == "and":
            children_converted = add_parens(children_elementary, children_converted)
            return elementary, "%s and %s " % (children_converted[0], children_converted[1])
        elif operator == "delay":
            if output_type in ["executable", "sympy"]:
                return elementary, children_converted[0]
            return elementary, "delay(%s, %s)" % (children_converted[0], children_converted[1])
        elif operator in ["exp", "ln", "floor", "ceiling", "factorial", "abs", "cos", "sin", "tan", "sinh", "cosh",
                          "tanh", "arcsin", "arccos", "arctan"]:
            return elementary, "%s(%s)" % (convert_function(output_type, operator), children_converted[0])
        elif operator == "log":
            if len(children_converted) == 1:
                return elementary, "%s(%s)" % (convert_function(output_type, "log"), children_converted[0])
            else:
                if output_type == "executable":
                    # NB. second argument to math.log is base [opposite to mathML]
                    return elementary, "math.log(%s, %s)" % (children_converted[1], children_converted[0])
                elif output_type == "sympy":
                    return elementary, "sympy.log(%s, %s)" % (children_converted[1], children_converted[0])

                return elementary, "%s_%s(%s)" % ("log", children_converted[0], children_converted[1])

        elif operator == "root":

            # default to sqrt()
            if len(children_converted) == 1:
                return elementary, "%s(%s)" % (convert_function(output_type, "sqrt"), children_converted[0])

            # otherwise root(n,a) = pow(a, 1/n)
            if len(children_converted) == 2:
                if output_type == "executable":
                    return elementary, "pow(%s, 1/%s)" % (children_converted[1], children_converted[0])
                return elementary, "%s(%s, %s)" % (convert_function(output_type, "root"), children_converted[0], children_converted[1])

    elif expression.name == "logbase":

        for child in expression.children:
            if isinstance(child, NavigableString):
                continue

            child_elementary, child_converted = convert_rate_law_inner(child, initial_values, non_default_variables, non_default_values, output_type)
            return child_elementary, child_converted

    if expression.name == "degree":
        # degree tag used with root
        for child in expression.children:
            if not isinstance(child, NavigableString):
                return convert_rate_law_inner(child, initial_values, non_default_variables, non_default_values, output_type)


def inline_all_functions(model):
    """
    Replace all uses of user-defined functions in a model with an inline use of the corresponding definition.

    This is a safe to perform: the SBML L3V1 Core Specification states:
    "With the restrictions as they are, function definitions could, if desired, be implemented as textual substitutions"

    Parameters
    ----------
    model : bs4.BeautifulSoup object produced by parsing an SBML model

    Returns
    -------
    model with uses of user-defined functions replaced
    """

    # Get function definitions
    function_definition = {}

    function_definitions = model.select_one('listOfFunctionDefinitions')

    if not function_definitions:
        return model

    for function in function_definitions.select("functionDefinition"):

        function_id = function.attrs["id"]
        math = copy.copy(function.select_one("math").select_one("lambda"))

        # get list of arguments to this function
        args = []
        for bvar in math.select('bvar'):
            args.append(bvar.select_one('ci').text.strip())

        # now remove the bvars the get the body of the function
        for bvar in math.select("bvar"):
            bvar.replace_with('')

        inner = ""
        for child in math.contents:
            if not isinstance(child, NavigableString):
                inner = child
                break

        function_definition[function_id] = {"math": inner, "arguments": args}

    # Replace function calls with inlined definitions
    replaced = True
    while replaced:
        replaced = False

        for math in model.select("math"):
            for apply_element in math.select('apply'):
                # get list of tag children
                children = []
                for child in apply_element.contents:
                    if not isinstance(child, NavigableString):
                        children.append(child)

                name = children[0].text.strip()
                if name in list(function_definition.keys()):
                    inlined = inline_function_call(function_definition[name], children[1:])
                    apply_element.replace_with(inlined)
                    replaced = True
                    break
    return model


def inline_function_call(func, arguments):
    """

    Parameters
    ----------
    func : dict representing user defined function
    arguments : BeautifulSoup object representing the expressions used as arguments to the function

    Returns
    -------
    BeautifulSoup object representing the supplied expressions substituted into the function definition
    """
    math = func["math"]
    args = func["arguments"]

    math = copy.copy(math)

    # for each arg, get list of ci elements
    cis = {}
    for ci in math.select("ci"):
        variable_name = ci.text.strip()
        if variable_name in args:
            if variable_name not in list(cis.keys()):
                cis[variable_name] = []
            cis[variable_name].append(ci)

    # now replace each element in this list
    for i in range(len(args)):
        if args[i] in list(cis.keys()):
            for ci in cis[args[i]]:
                ci.replace_with(copy.copy(arguments[i]))

    return math
