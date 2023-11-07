import numpy as np

def lu_decomposition_with_row_pivoting(A):
    n = len(A)
    P = np.eye(n)
    for i in range(n):
        p = np.argmax(abs(A[i:, i])) + i
        if p != i:
            A[[i, p]] = A[[p, i]]
            P[[i, p]] = P[[p, i]]
        for j in range(i+1, n):
            A[j, i] /= A[i, i]
            A[j, i+1:] -= A[j, i] * A[i, i+1:]
    L = np.tril(A, k=-1) + np.eye(n)
    U = np.triu(A)
    return L, U, P

def solve_system_with_lu_decomposition(A, b):
    L, U, P = lu_decomposition_with_row_pivoting(A)
    y = np.linalg.solve(L @ P @ A, P @ b)
    x = np.linalg.solve(U @ P @ A, y)
    return x

def inverse_matrix(A):
    n = len(A)
    I = np.eye(n)
    X = np.zeros((n, n))
    for i in range(n):
        X[:, i] = solve_system_with_lu_decomposition(A, I[:, i])
    return X

# Пример нахождения обратной матрицы для матрицы Лемера размерности n=3
A = np.array([[1, 1/2, 1/3], [1/2, 2/3, 1], [1/3, 1, 3/4]])
A_inv = inverse_matrix(A)
print('Обратная матрица для матрицы Лемера (n = 3):\n', A_inv)
print("A * A^-1:\n", np.round(A @ A_inv, decimals=5))
