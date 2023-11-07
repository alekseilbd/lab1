import numpy as np
import matplotlib.pyplot as plt

# Краевые условия
T_left = 0
T_right = 100
L = 1  # Длина стержня
t_max = 1  # Время

# Шаг сетки
dx = 0.01
dt = 0.001

# Число узлов
N = int(L / dx) + 1
M = int(t_max / dt) + 1

# Константы
alpha = 1  # Коэффициент теплопроводности
r = alpha * dt / dx**2

# Инициализация сетки
T = np.zeros((M, N))
T[0, :] = T_left
T[-1, :] = T_right

# Цикл по времени
for m in range(1, M):
    # Метод переменных направлений
    A = np.eye(N) * (1 + 2 * r)
    A[0, 1] = -r
    A[-1, -2] = -r
    B = np.eye(N) * (1 - 2 * r)
    B[0, 1] = r
    B[-1, -2] = r
    T[m, :] = np.linalg.solve(A, np.dot(B, T[m-1, :]))

# График решения
x = np.linspace(0, L, N)
t = np.linspace(0, t_max, M)
X, T = np.meshgrid(x, t)
fig = plt.figure(figsize=(8, 6))
ax = fig.gca(projection='3d')
ax.plot_surface(X, T, T, cmap='coolwarm')
ax.set_xlabel('X')
ax.set_ylabel('Time')
ax.set_zlabel('Temperature')
plt.show()
