import numpy as np
from scipy.linalg import lu

n = int(input("Введите размерность матрицы: "))
A = [[min(g, j) for j in range(1, n+1)] for g in range(1, n+1)]

P, L, U = lu(A)

A_inv = np.linalg.inv(U) @ np.linalg.inv(L) @ np.linalg.inv(P)

norm_A = np.linalg.norm(A, ord=np.inf)
norm_A_inv = np.linalg.norm(A_inv, ord=np.inf)

cond_A = norm_A * norm_A_inv
print("Число обусловленности матрицы А:", cond_A)
