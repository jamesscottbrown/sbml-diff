from bs4 import NavigableString
from rate_laws import convert_rate_law
import math  # needed for check_sign_numerically()


def categorise_interaction(kinetic_law, species_id):
    for math_expr in kinetic_law.select_one("math"):
        if isinstance(math_expr, NavigableString):
            continue
        return categorise_interaction_inner_numerically(math_expr, species_id)


def categorise_interaction_inner_numerically(expression, species_id):
    # identify all parameters and concentrations in the rate law
    symbols = []
    for ci in expression.findAll("ci"):
        symbols.append(ci.text.strip())
    symbols = set(symbols)

    return check_sign_numerically(expression, symbols, species_id)


def check_sign_numerically(expr, param_names, species_id):
    # This function will fail to report that a kineticLaw is mixed (rather than monotonic) if:
    # - the sign of its gradient depends on value of one of the parameters (eg x^a/x^b)
    # - the sign of its gradient depends on the concentration of the corresponding species

    expr = convert_rate_law(expr, variables_not_to_substitute=[species_id])

    replacement = {'^': '**', 'exp': 'math.exp', 'ln': 'math.log', 'log': 'math.log10', 'ceiling': 'math.ceil',
                   'floor': 'math.floor', 'factorial': 'math.factorial', 'delay': '', 'pi': 'math.pi',
                   'exponentiale': 'math.e', 'infinity': 'float("Inf")'}

    for old in replacement:
        expr = expr.replace(old, replacement[old])

    param_names.remove(species_id)

    rate_change = eval(expr.replace(species_id, '1')) - eval(expr.replace(species_id, '0'))

    if rate_change > 0:
        return "monotonic_increasing"
    elif rate_change == 0:
        return "constant"
    else:
        return "monotonic_decreasing"
