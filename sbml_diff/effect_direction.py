from bs4 import NavigableString, Tag
from .rate_laws import convert_rate_law
import math  # needed for check_sign_numerically()


def categorise_interaction(kinetic_law, species_id, initial_values, use_sympy=False):
    """
    Given a kineticLaw and the name of a species, determine whether the expression is a monotonic_increasing,
    monotonic_decreasing, or constant with respect to the concentration of that species.

    Parameters
    ----------
    kinetic_law : bs4.element.Tag corresponding to a kineticLaw
        
    species_id : the species id
        

    Returns
    -------
    string representing the sign of the interaction

    """
    for math_expr in kinetic_law.select_one("math"):
        if isinstance(math_expr, NavigableString):
            continue

        # identify all parameters and concentrations in the rate law
        symbols = []

        if isinstance(math_expr, Tag) and math_expr.name == "ci":
            symbols.append(math_expr.text.strip())
        else:
            for ci in math_expr.findAll("ci"):
                symbols.append(ci.text.strip())
        symbols = set(symbols)

        if use_sympy:
            return check_sign_algebraically(math_expr, symbols, species_id, initial_values)
        else:
            return check_sign_numerically(math_expr, symbols, species_id, initial_values)


def check_sign_algebraically(expr, param_names, species_id, initial_values):

    import sympy
    from sympy import Symbol

    # create variables for each param and species
    for parameter_name in param_names:
        exec("%s = Symbol('%s', positive=True)" % (parameter_name, parameter_name))

    # create variables for time and avogadro's constant: N_A
    exec("t = Symbol('t', positive=True)")
    exec("N_A = Symbol('N_A', positive=True)")

    exec("query_var = %s" % species_id)

    converted_rate_law = convert_rate_law(expr, output_type="sympy")
    if converted_rate_law == "piecewise":
        return "?"

    exec("rate = %s" % converted_rate_law)
    if sympy.simplify(rate.diff(query_var)).is_positive:
        return "monotonic_increasing"
    elif sympy.simplify(rate.diff(query_var)).is_negative:
        return "monotonic_decreasing"
    else:
        return "?"

def check_sign_numerically(expr, param_names, species_id, initial_values):
    """
    Given a MathML expression, list of all parameter/species names, and the name of a species, determine whether the
    expression is a monotonic_increasing, monotonic_decreasing, or constant with respect to the concentration of that
    species.

    This is done by setting all parameters and the concentrations of other species to 1, then comparing the value of the
    expression when the concentration of interest is 1 or 0.
    This approach will fail to report that a kineticLaw is mixed (rather than monotonic) if:
    - the sign of its gradient depends on value of one of the parameters (eg x^a/x^b)
    - the sign of its gradient depends on the concentration of the corresponding species


    Parameters
    ----------
    expr : bs4.element.Tag object corresponding to the contents of a math element
        
    param_names : list of the names of all parameters
        
    species_id : list of the species we are interested in
        

    Returns
    -------

    """

    expr1 = convert_rate_law(expr, initial_values, non_default_variables=[species_id], non_default_values='1', output_type="executable")
    expr2 = convert_rate_law(expr, initial_values, non_default_variables=[species_id], non_default_values='0.01', output_type="executable")

    if not expr1 or not expr2:
        return '?'

    param_names.remove(species_id)

    try:
        rate_change = eval(expr1) - eval(expr2)
    except (ZeroDivisionError, ValueError, OverflowError):
        return "?"

    if rate_change > 0:
        return "monotonic_increasing"
    elif rate_change == 0:
        return "constant"
    else:
        return "monotonic_decreasing"
