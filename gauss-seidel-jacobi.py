import math
import numpy as np
import pandas as pd

dat = {}
dat2 = {}
errors = {}
errors2 = {}

def diagDom(A, n):
    for i in range(n):
        max = A[i, i]
        sum = 0
        for j in range(n):
            if i != j:
                sum += A[j][i]
        if sum > max:
            return -1
    return 1


def jacobi(A, b, x0, e, N):
    n = np.shape(x0)[0]
    x1 = np.zeros(n)
    k = 0
    dat[0] = np.copy(x0)
    errors[0] = [0.0]
    while k < N:
        # print(x0)
        k += 1
        for i in range(0, n):
            sum = 0
            for j in range(0, n):
                if j != i:
                    sum += (A[i][j] * x0[j])
            x1[i] = (b[i] - sum) / A[i][i]

        error = np.linalg.norm(x1 - x0)
        dat[k] = np.copy(x1)
        errors[k] = [error]
        x0 = np.copy(x1)
        if error < e:
            # print(error)
            return x1, k

    if k > N:
        return "Not found"


def gauss_seidel(A, b, x, e, N):
    k = 0
    dat2[0] = np.copy(x)
    errors2[0] = [0.0]
    while k < N:
        k += 1
        s = 0.0

        for i in range(0, n):
            t = x[i]
            sum = 0
            for j in range(0, n):
                if i != j:
                    sum += (A[i][j] * x[j])
            x[i] = (b[i] - sum) / A[i][i]
            s += ((x[i] - t) ** 2)

        error = math.sqrt(s)
        errors2[k] = [error]
        dat2[k] = np.copy(x)
        if error < e:
            # print(error)
            return x, k

    return "Not found"


if __name__ == "__main__":
    pd.set_option("display.max_rows", None, "display.max_columns", None, 'display.width', 10000000)

    A = np.array(np.mat('8 4 -3; 1 -10 2; 3 3 -7'))
    n = np.shape(A)[0]
    x0 = np.zeros(3)
    b = np.array(np.mat('3;2;1'))
    # A = np.array(np.mat('10 -2 8; 3 6 4; 5 3 13'))
    # n = np.shape(A)[0]
    #
    # b = np.array(np.mat('1;2;3'))
    # x0 = np.ones(3)
    print(A)
    print(b)

    # print(diagDom(A, n))
    e = 10 ** -5
    print("Jacobi\n")
    ans = jacobi(A, b, x0, e, 50)
    e1 = pd.DataFrame.from_dict(errors)
    df = pd.DataFrame.from_dict(dat).append(e1)

    print(df)
    print('ans = {0}\nIterations = {1}'.format(ans[0],ans[1]))
    print("\n")
    print("GAUSS SEIDEL\n")
    ans = gauss_seidel(A, b, x0, e, 25)
    e2 = pd.DataFrame.from_dict(errors2)
    df2 = pd.DataFrame.from_dict(dat2).append(e2)
    print(df2)
    print('ans = {0}\nIterations = {1}'.format(ans[0],ans[1]))
    # print(A)
    # print(b)
    # print(x0)
