from sympy import sin, cos, Function, Symbol, exp, Pow
from sympy.abc import l, theta, phi, beta, i, m

from print_helpers import print_equation


def hyperspharm():
    Znlm = Function("Z__n_l,m")
    NPml = Function("NP__m_l")
    NGml = Function("NG__l+1_n-l")
    print_equation(Znlm(beta, theta, phi), NPml(cos(theta)) * NGml(cos(beta)) * Pow(sin(theta), l) * exp(i * m * phi))


if __name__ == '__main__':
    hyperspharm()
