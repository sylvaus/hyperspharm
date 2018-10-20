from sympy import sin, cos, Function, exp, Pow, Integral, pi
from sympy.abc import l, theta, phi, beta, i, m, x

from print_helpers import print_equation
import legendre as lg
import gegenbauer as gg


def NPml(m_, l_, Pmlx):
    return lg.spharm_norm_coeff(m_, l_) * Pmlx


def NGml(m_, l_, Gmlx):
    return gg.norm_coeff(m_, l_) * Gmlx


def integral_over_hyper_sphere(func):
    return Integral(Integral(Integral(func, (beta, 0, pi)),(theta, 0, pi)), (phi, 0, 2*pi))

def hyperspharm():
    Znlm_name = Function("Z__m_l,n")
    NPml_name = Function("NP__m_l")
    NG_l_1_n_l_name = Function("NG__l+1_n-l")
    print_equation(Znlm_name(beta, theta, phi),
                   NPml_name(cos(theta)) * NG_l_1_n_l_name(cos(beta)) * Pow(sin(theta), l) * exp(i * m * phi))

    print_equation(NPml_name(x), NPml(m, l, Function("P__m_l")(x)), True)
    print_equation(Function("NG__m_l")(x), NGml(m, l, Function("G__m_l")(x)), True)


def hyperspharm_transform():
    Znlm_name = Function("Z__m*_l,n")
    fnlm = Function("f__m_l,n")
    f_func = Function("f")
    print_equation(
        fnlm,
        integral_over_hyper_sphere(f_func(beta, theta, phi) *
                                   Znlm_name(beta, theta, phi) *
                                   sin(beta) ** 2 * sin(theta)),
        True
    )


if __name__ == '__main__':
    hyperspharm()
    hyperspharm_transform()
