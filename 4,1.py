

import numpy as np


def LU_decomposition(A):
    n = len(A)
    L = np.zeros((n, n))
    U = np.zeros((n, n))
    P = np.eye(n)

    for i in range(n):
        max_row = i
        for j in range(i + 1, n):
            if abs(A[j, i]) > abs(A[max_row, i]):
                max_row = j
        if max_row != i:
            A[[i, max_row]] = A[[max_row, i]]
            P[[i, max_row]] = P[[max_row, i]]

        L[i, i] = 1
        for j in range(i):
            L[i, j] = A[i, j]
            for k in range(j):
                L[i, j] -= L[i, k] * U[k, j]
            L[i, j] /= U[j, j]

        for j in range(i, n):
            U[i, j] = A[i, j]
            for k in range(i):
                U[i, j] -= L[i, k] * U[k, j]

    return L, U, P


def dec_LU(A):
     n=len(A)
     LU=np.copy(A)
     for j in range(0,n-1):
         for i in range(j+1,n):
             if LU[i,j]!=0:
                 u=(LU[i,j])/LU[j,j]
                 LU[i,j+1:n]=LU[i,j+1:n]-u*LU[j,j+1:n]
                 LU[i,j]=u
     return LU

n = 4
Pascal = np.zeros((n, n))
for i in range(n):
    Pascal[i, 0] = 1
    Pascal[0, i] = 1
for i in range(1, n):
    for j in range(1, n):
        Pascal[i, j] = Pascal[i - 1, j] + Pascal[i, j - 1]
print("Matrix Pascal:")
print(Pascal)

L, U, P = LU_decomposition(Pascal)
LU_=np.zeros((n,n))
for i in range(0, n):
    for j in range(0, n):
        for k in range(0, n):
            LU_[i][j] += L[i][k] * U[k][j]
print("Матрица L:")#L — нижняя треугольная матрица с единичной диагональю
print(L)
print("МАтрица U:")#— верхняя треугольная матрица с ненулевыми диагональными элементами

print(U)
print("Матрица P:")#Матрица перестановок
print(P)

print("LU разложение модуль lu:")
print(dec_LU(Pascal))
