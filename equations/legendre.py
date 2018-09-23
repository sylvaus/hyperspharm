from sympy import pi, factorial, sqrt

def spharm_norm_coeff_square(m_, l_):
    return (1 / (4 * pi)) * ((2 * l_ + 1) * (factorial(l_ - m_)) / factorial(l_ + m_))


def spharm_norm_coeff(m_, l_):
    return 1 / sqrt(4 * pi) * sqrt((2 * l_ + 1) * (factorial(l_ - m_)) / factorial(l_ + m_))
