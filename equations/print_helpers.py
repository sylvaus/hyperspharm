from sympy import pprint, Equality, simplify, latex


def print_equation(before_equal, after_equal, no_simplify=False):
    if not no_simplify:
        after_equal = simplify(after_equal)
        before_equal = simplify(before_equal)

    pprint(Equality(before_equal, after_equal))
    print(latex(Equality(before_equal, after_equal)))
