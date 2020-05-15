import math
import numpy as np
import pandas as pd
from numpy.polynomial.polynomial import Polynomial
from scipy.interpolate import lagrange


def getLagrange(x, y, x0=None, roundTo=None):
    print("\nLAGRANGE\n")
    print("X = " + str(x))
    print("Y = " + str(y))
    poly = Polynomial(lagrange(x, y)).coef
    print(poly)
    if roundTo:
        poly = list(map(lambda x: round(x, roundTo), poly))
    print(coefToNiceString(poly[::-1]))

    if x0:
        n = len(poly) - 1
        sum = poly[n]
        for i in reversed(range(0, n)):
            # print('{0} {1}'.format(n - i, poly[i]))
            sum += poly[i] * (x0 ** (n - i))
        print("\n Approximation using x0: {} = {}".format(x0, sum))


def coefToNiceString(coef):
    coef = coef[::-1]
    res = ""
    degree = len(coef) - 1
    res += str(coef[0]) + "x^" + str(degree)
    for i in range(1, len(coef) - 1):
        coeff = coef[i]
        if coeff < 0:
            res += " - " + str(-coeff) + "x^" + str(degree - i)
        else:
            res += " + " + str(coeff) + "x^" + str(degree - i)

    if coef[-1] < 0:
        res += " - " + str(-coef[-1])
    else:
        res += " + " + str(coef[-1])

    return res


def create_matrix(x, y):
    n = len(y) if len(y) >= len(x) else len(x)
    return np.eye(n) * y, n


def compute_divided_differences(x, y, x0=None, roundTo=None):
    print("\nDIVIDED DIFFERENCES")
    print("X = " + str(x))
    print("Y = " + str(list(y)))
    print("\n")
    A, n = create_matrix(x, y)

    for d in range(1, n):
        for i in range(0, n - d):
            j = i + d
            A[i][j] = (A[i + 1][j] - A[i][j - 1]) / (x[j] - x[i])
            print('f[{0}][{1}] = ({2} - {3})/({4} - {5})'.format(i, j, A[i + 1][j], A[i][j - 1], x[j], x[i]))
    print("\n")
    print("X = " + str(x))
    print(A)
    print("\n")
    print("=> Polynomial: {}".format(dividedDiffPoly(A[0], x)))
    if x0:
        print('Poly: degree {0} | x0 = {1} | approximation = {2}'.format(n - 1, x0, approx_div_diff(A[0], x, x0)))

    return A[0]


def getSign(x):
    return " + " if x > 0 else " - "


def dividedDiffPoly(coefs, x):
    poly = "{0}".format(coefs[0])
    for i in range(1, len(coefs)):
        calcXs = ""
        for j in range(i):
            sign = getSign(-1 * x[j])
            calcXs += "(x{0}{1})".format(sign, abs(x[j]))
        poly += ' {0} {1}*{2}'.format(getSign(coefs[i]), abs(coefs[i]), calcXs)
    # print(poly)
    return poly


def approx_div_diff(coefs, x, x0):
    sum = coefs[0]
    for i in range(1, len(coefs)):
        calcXs = 1
        for j in range(i):
            calcXs *= (x0 - x[j])
        sum += (coefs[i] * calcXs)
    return sum


def compute_nevilles(x, y, x0, roundTo=None):
    print("\nNEVILLES\n")
    print("X = " + str(x))
    print("Y = " + str(y))
    A, n = create_matrix(x, y)

    for d in range(1, n):
        for i in range(0, n - d):
            j = i + d
            A[i][j] = (((x0 - x[i]) * A[i + 1][j]) - ((x0 - x[j]) * A[i][j - 1])) / (x[j] - x[i])

            print('P[{0}][{1}] = (({2}-{3})*{4}) - ({2}-{5})*{6}))/({5}-{3}) => {7}'.format(i, j, x0, x[i], A[i + 1][j],
                                                                                            x[j],
                                                                                            A[i][j - 1], A[i][j]))
    print("\n")
    print("X = " + str(x))
    # print("Y = " + str(y))
    print(A)
    print("\n Approximation using x0: {} = {}".format(x0, A[0, n - 1]))

    return A[0, n - 1]


if __name__ == "__main__":
    # x = [1, 1.3, 1.6, 1.9, 2.2]
    # y = [0.7651977, 0.6200860, 0.4554022, 0.2818186, 0.1103623]
	
    # x = [1, 2, 3, 4, 5]
    # y = [.7, .73, .8, .75, .6]
	
    # x = [1,1.3,1.6,1.9]
    # y = [0.7651977,0.6200860,0.4554022,0.2711463]
	
    # x = [0, 1, 2, 3]
    # y = [2.8, 3.5, 1.6, 3.0]
	
    # x = [1.1, 1.2, 1.3, 1.4]
    # y = [2.3, 2.4, 2.1, 3.6]
	
    # x = [-1, 1, 3, 7]
    # y = [18, -8, -30, 130]
    # x = [1, 2, 3, 4]
    # y = np.zeros(len(x))
    # for i in range(len(y)):
        # y[i] = x[i] ** 4 + math.sqrt(2) * (x[i] ** 3) + math.pi * x[i]
		
    # x = [-25, -5, 5, 10, 20, 30]
    # y = [424585, 1025, 445, 8330, 145450, 759370]
	
    x = [0.6283185308, 1.2566370616, 1.8849555924, 2.5132741232]
    y = [0.587785252358846, 0.859485618518891, 0.951056516219097, 0.587785252026982]
    x0 = 1.5
    

    # DIVIDED DIFFERENCES #


    # print("Can be reduced to: 11.41421356237310x^3 - 35x^2 + 53.1415926535899x - 24")
    #
    # x = [0, 1, math.pi, math.e, -1]
    # y = np.zeros(len(x))
    # for i in range(len(y)):
    #     y[i] = x[i] ** 4 + math.sqrt(2) * (x[i] ** 3) + math.pi * x[i]
    # sol = compute_divided_differences(x, y)

    print("\n---The polynomials are different degrees so they are not the same---")
    # print("Can be reduced to: x^4 + 1.41421x^3 + 7.10543×10^-15x^2 + 3.14159x")
    # NEVILLES #

    x = [0,1,2,3]
    y = [2.8,3.5,1.6,3.0]

    compute_nevilles(x, y, 3.5, 0)
    compute_divided_differences(x, y,3.5)
    # x0 = .9
    # x = [.5, .6, .7, .8, 1]
    # y = [-0.34409873, -0.17694460, 0.01375227, 0.22363362, 0.65809197]

    getLagrange(x, y, 3.5, 8)
