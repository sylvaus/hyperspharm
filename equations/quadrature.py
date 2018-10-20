from typing import List
from math import cos, pi
from collections import namedtuple

FunctionResult = namedtuple("FunctionResult", ["func", "result", "str"])

FUNCTION_RESULT_N_4 = [
    FunctionResult(lambda x: x, 0, "f(x) = x")
    , FunctionResult(lambda x: x ** 2, 2/3, "f(x) = x ** 2")
    , FunctionResult(lambda x: x ** 3, 0, "f(x) = x ** 3")
    , FunctionResult(lambda x: x ** 4, 2/5, "f(x) = x ** 4")
]

FUNCTION_RESULT_N_13 = [
    FunctionResult(lambda x: x, 0, "f(x) = x")
    , FunctionResult(lambda x: x ** 10, 2/11, "f(x) = x ** 10")
    , FunctionResult(lambda x: x ** 11, 0, "f(x) = x ** 11")
    , FunctionResult(lambda x: x ** 12, 2/13, "f(x) = x ** 12")
    , FunctionResult(lambda x: (x ** 12) +  2 *(x ** 10) + 3 * (x ** 8) + (x**5) + 5 * x,
                     2 * (1/13 + 2/11 + 3/9),
                     "f(x) = (x ** 12) +  2 *(x ** 10) + 3 * (x ** 8) + (x**5) + 5 * x")
]


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def trapezoid_rule(points: List[Point]):
    result = 0
    for i in range(len(points) - 1):
        result += (points[i + 1].y + points[i].y) * abs(points[i + 1].x - points[i].x) / 2.0

    return result


def clenshaw_curtis_rule_even_only(points: List[Point]):
    """
    Based on https://en.wikipedia.org/wiki/Clenshaw%E2%80%93Curtis_quadrature
    """
    N = (len(points) - 1)
    d = [1.0] + \
        [2.0 / (1.0 - (2.0 * k) ** 2) for k in range(1, N // 2)] + \
        [1.0 / (1.0 - N ** 2)]
    result = 0
    ws = []
    for n in range((N // 2) + 1):
        w = 0
        for k, dj in enumerate(d):
            if (n == 0) or (n == (N // 2)):
                mult = 0.5
            else:
                mult = 1.0
            w += mult * (2.0 / N) * cos(n * k * pi * 2.0 / N) * dj
        ws.append(w)
        result += ((points[n].y + points[-(1 + n)].y) * w)

    return result


def clenshaw_curtis_rule_expanded_even_only(points: List[Point]):
    N = (len(points) - 1)
    d = [1.0] + \
        [2.0 / (1.0 - (2.0 * k) ** 2) for k in range(1, N // 2)] + \
        [1.0 / (1.0 - N ** 2)]
    result = 0
    ws = []
    for n in range(N + 1):
        w = 0
        for k, dj in enumerate(d):
            if (n == 0) or (n == N):
                mult = 0.5
            else:
                mult = 1.0
            w += mult * (2.0 / N) * cos(n * k * pi * 2.0 / N) * dj
        ws.append(w)
        result += (points[n].y * w)

    return result

def clenshaw_curtis_rule_expanded(points: List[Point]):
    N = (len(points) - 1)
    if (len(points) % 2) == 0:
        N2max = (N - 1) // 2
    else:
        N2max = N // 2

    result = 0
    ws = []
    for n, point in enumerate(points):
        w = 0
        for k in range(N2max + 1):
            # c2k
            if (k == 0) or (2*k == N):
                c2k = 1.0 / (1.0 - (2.0 * k) ** 2)
            else:
                c2k = 2.0 / (1.0 - (2.0 * k) ** 2)
            # bn
            if (n == 0) or (n == N):
                bn = 0.5
            else:
                bn = 1.0
            w += bn * (2.0 / N) * cos(n * k * pi * 2.0 / N) * c2k
        ws.append(w)
        result += (point.y * w)

    return result

def compare_str(value, expected):
    return "value: {}, abs error: {}".format(value, abs(value - expected))

def verification(n, function_results):
    for func_result in function_results:
        points = []
        for i in range(n + 1):
            x = cos(pi * (i / n))
            points.append(Point(x, func_result.func(x)))

        print("\nTest with {}."
              "\nExpected result:                   {}".format(func_result.str, func_result.result))
        print("Trapezoid                         ",
              compare_str(trapezoid_rule(points), func_result.result))
        print("Clenshaw Curtis even only factored",
              compare_str(clenshaw_curtis_rule_even_only(points), func_result.result))
        print("Clenshaw Curtis even only expanded",
              compare_str(clenshaw_curtis_rule_expanded_even_only(points), func_result.result))
        print("Clenshaw Curtis expanded          ",
              compare_str(clenshaw_curtis_rule_expanded(points), func_result.result))


if __name__ == '__main__':
    verification(13, FUNCTION_RESULT_N_13)
