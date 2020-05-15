import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

dat = {'xi': [], 'hi': [], 'yi=ai': [], 'ci': [], 'bi': [], 'di': []}


def solve_for_c(x, y):
    n = len(x)
    h = np.zeros(n - 1)
    b = np.zeros(n)
    A = np.eye(n)

    for i in range(n - 1):
        h[i] = (x[i + 1] - x[i])
    # print("h = " + str(h))
    dat['hi'] = h

    for i in range(1, n - 1):
        A[i, i - 1] = h[i - 1]
        A[i, i] = 2 * (h[i - 1] + h[i])
        A[i, i + 1] = h[i]
        # breakpoint()
        b[i] = 3 * (((y[i + 1] - y[i]) / h[i]) - ((y[i] - y[i - 1]) / h[i - 1]))
    print("Matrix of h values = \n" + str(A))
    b = np.reshape(b, (n, 1))
    print("b = " + str(b))
    c = np.matmul(np.linalg.inv(A), b)
    return c.flatten()


def solve_for_d(x, c):
    n = len(x)
    d = np.zeros(n - 1)
    for i in range(n - 1):
        h = (x[i + 1] - x[i])
        d[i] = ((c[i + 1] - c[i]) / (3 * h))
    return d


def solve_for_b(x, y, c):
    n = len(x)
    b = np.zeros(n - 1)
    for i in range(n - 1):
        h = (x[i + 1] - x[i])
        b[i] = ((y[i + 1] - y[i]) / h) - h * (2 * c[i] + c[i + 1]) / 3
    return b


def make_func(start, end, a, b, c, d):
    x = np.arange(start, end + 0.05, 0.05)
    y = a + b * (x - start) + c * (x - start) ** 2 + d * (x - start) ** 3
    return x, y


def plot_splines(x, y, c, d, b):
    fig = plt.figure()
    plt.plot(x, y, 'o')
    n = len(x)
    for i in range(n - 1):
        xi, yi = make_func(x[i], x[i + 1], y[i], b[i], c[i], d[i])
        plt.plot(xi, yi)
    plt.show()


if __name__ == "__main__":

    print("\n#5) ")
    x = [0.025, 0.4, 0.6, .75, 1]
    #y = np.zeros(len(x))
    x = [1, 2, 3, 4, 5, 6]
    y = [5, 6, 6.5, 5.5, 5.5, 7]
     #x = list(range(1, 7))
    # y = [5, 6, 6.5, 5.5, 5.5, 7]
    # print(x)
    # print(y)

    c = solve_for_c(x, y)
    d = solve_for_d(x, c)
    b = solve_for_b(x, y, c)
    plot_splines(x, y, c, d, b)

    dat['xi'] = x
    dat['yi=ai'] = y
    dat['ci'] = c
    dat['bi'] = b
    dat['di'] = d

    df = pd.DataFrame.from_dict(dat, orient='index')
    df = df.transpose()
    print(df)

    for index, row in df.iterrows():
        if index != len(x)-1:
            print('S{0} = {1} + {2}(x-{3}) + {4}(x-{3})^2 + {5}(x-{3})^3'.format(index, row['yi=ai'], row['bi'],
                                                                                 row['xi'],
                                                                                 row['ci'], row['di']))
