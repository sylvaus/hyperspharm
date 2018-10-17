from typing import List
from math import cos, pi


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y


def trapezoid_rule(points: List[Point]):
    result = 0
    for i in range(len(points) - 1):
        result += (points[i + 1].y + points[i].y) * abs(points[i + 1].x - points[i].x) / 2.0

    return result


def clenshaw_curtis_rule(points: List[Point]):
    N = (len(points) - 1)
    d = [1.0] + \
        [2.0 / (1.0 - (2.0 * k) ** 2) for k in range(1, (N // 2) - 1)] + \
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
        result += ((points[n].y + points[-(1+n)].y) * w)

    return result


def verification():
    def func(x):
        return x ** 4 + x ** 3 + x ** 11

    N = 12
    points = []
    for i in range(N + 1):
        x = cos(pi * (i / N))
        points.append(Point(x, func(x)))

    print("Trapezoid", trapezoid_rule(points))
    print("Clenshaw Curtis", clenshaw_curtis_rule(points))


if __name__ == '__main__':
    verification()
