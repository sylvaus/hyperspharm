from sympy import Function, simplify, sqrt, factorial, pi, Symbol
from sympy.abc import x, m, l

from print_helpers import print_equation


def aml(m_, l_):
    return (2 * x * (m_ + l_ - 1)) / l_


def bml(m_, l_):
    return - (2 * m_ + l_ - 2) / l_


def norm_coeff_square(m_, l_):
    return ((2 ** (2 * m_ - 1) * (m_ + l_) * factorial(l_))
            / (pi * factorial(2 * m_ + l_ - 1))) \
           * factorial(m_ - 1) ** 2


def gegenbauer():
    Gml = Function("G__m_l")
    Gml_1 = Function("G__m_l-1")
    Gml_2 = Function("G__m_l-2")

    Gm0 = Function("G__m_0")
    Gm1 = Function("G__m_1")

    print_equation(Gm0, Symbol("1"))
    print_equation(Gm1, 2 * m * x)

    print_equation(Gml(x), simplify(aml(m, l) * Gml_1(x) + bml(m, l) * Gml_2(x)))

    print_equation(Function("N__1_0")(x), sqrt(simplify(norm_coeff_square(1, 0))))

    print_equation(Function("N__1_1")(x), 2 * x * sqrt(simplify(norm_coeff_square(1, 1))))

    print_equation((Function("Factor__m+1_0")(x) / Function("Factor__m_0")(x)) ** 2,
                   norm_coeff_square(m + 1, 0) / norm_coeff_square(m, 0))

    print_equation((Function("Factor__m+1_1")(x) / Function("Factor__m_1")(x)) ** 2,
                   norm_coeff_square(m + 1, 1) / norm_coeff_square(m, 1))

    print_equation(Function("N__m+1_0")(x),
                   sqrt(simplify(norm_coeff_square(m + 1, 0) / norm_coeff_square(m, 0))) * Function("N__m_0")(x))

    print_equation(Function("N__m+1_1")(x),
                   sqrt(simplify(
                       (norm_coeff_square(m + 1, 1) / norm_coeff_square(m, 1)) * (((m + 1) / m) ** 2))) * Function(
                       "N__m_1")(x))

    print_equation(Function("N__m_l")(x),
                   (sqrt(simplify(norm_coeff_square(m, l) / norm_coeff_square(m, l - 1) * (aml(m, l) ** 2)))
                    * Function("N__m_l-1")(x))
                   -  # The - sign comes from the negative bml being squared
                   (sqrt(simplify(norm_coeff_square(m, l) / norm_coeff_square(m, l - 2) * (bml(m, l) ** 2)))
                    * Function("N__m_l-2")(x)))


if __name__ == '__main__':
    gegenbauer()
