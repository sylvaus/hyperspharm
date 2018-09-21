from sympy import pprint, Equality, simplify, latex


def print_equation(before_equal, after_equal):
    pprint(Equality(before_equal, simplify(after_equal)))
    print(latex(Equality(before_equal, simplify(after_equal))))