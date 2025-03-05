from numpy import zeros, linspace, log, sqrt, sum, exp
import numpy as np
import matplotlib.pyplot as plt

# x^2y''  - 4xy' + (6-x^2)y = 0

# Функция f подготавливает массив, содержащий элементы вектор-функции,

def f(u,x):
    f = zeros(2)
    f[0] = u[1]
    f[1] = u[1]*4/x - u[0]*(6-x*x)/(x*x)
    return f

M = 20
M_r = 1000
t_0 = 1
T = 10
tau = (T - t_0)/M
tau_r = (T- t_0)/M_r
t = linspace(t_0, T, M+1)
tm = linspace(t_0, T, M_r+1)
u = zeros((M+1, 2))
u_r = zeros((M_r+1, 2))
y_0 = 1
z_0 = -1
u[0] = [y_0, z_0]
u_r[0] = [y_0, z_0]
        
for m in range(M):
    u[m+1] = u[m] + tau*f(u[m], t[m])

for m in range(M_r):
    u_r[m+1] = u_r[m] + tau_r*f(u_r[m], tm[m])


y = [-exp(x-1)*x*x+2*exp(-x+1)*x*x for x in tm]

plt.plot(tm, y, color='green', label='Аналитическое решение')
plt.plot(t, u[:, 0], color='blue', label='Эйлер, M=20')
plt.plot(tm, u_r[:, 0], color='red', label='Эйлер, M = 1000')
plt.xlim(1, 10)
plt.ylim(-10, 5)
plt.legend()
plt.title('Схема Эйлера')
plt.show()

u[0] = [y_0, z_0]
u_r[0] = [y_0, z_0]

for m in range(M):
    w_1 = f(u[m], t[m])
    w_2 = f(u[m] +tau*2/3*w_1, t[m] +tau*2/3)
    u[m+1] = u[m] + tau*(1/4*w_1 + 3/4*w_2)

for m in range(M_r):
    w_1 = f(u_r[m], tm[m])
    w_2 = f(u_r[m] +tau_r*2/3*w_1, tm[m] +tau_r*2/3)
    u_r[m+1] = u_r[m] + tau_r*(1/4*w_1 + 3/4*w_2)


plt.plot(t, y, color='green', label='Аналитическое решение')
plt.plot(t, u[:, 0], color='blue', label='ERK2, M=20')
plt.plot(tm, u_r[:, 0], color='red', label='ERK2, M = 1000')
plt.xlim(1, 10)
plt.ylim(-10, 5)
plt.legend()
plt.title('ERK2')
plt.show()

u[0] = [y_0, z_0]
u_r[0] = [y_0, z_0]
for m in range(M):
    w_1 = f(u[m], t[m])
    w_2 = f(u[m] + tau*1/2*w_1, t[m] + tau*1/2)
    w_3 = f(u[m] + tau*1/2*w_2, t[m] + tau*1/2)
    w_4 = f(u[m] + tau*1*w_3, t[m] + tau*1)
    u[m+1] = u[m] + tau*(1/6*w_1 + 1/3*w_2 + 1/3*w_3 + 1/6*w_4)

for m in range(M_r):
    w_1 = f(u_r[m], tm[m])
    w_2 = f(u_r[m] + tau_r*1/2*w_1, tm[m] + tau_r*1/2)
    w_3 = f(u_r[m] + tau_r*1/2*w_2, tm[m] + tau_r*1/2)
    w_4 = f(u_r[m] + tau_r*1*w_3, tm[m] + tau_r*1)
    u_r[m+1] = u_r[m] + tau_r*(1/6*w_1 + 1/3*w_2 + 1/3*w_3 + 1/6*w_4)

plt.plot(t, y, color='green', label='Аналитическое решение')
plt.plot(t, u[:, 0], color="blue", label='ERK4, M=20')
plt.plot(tm, u_r[:, 0], color='red', label='ERK4, M=1000')
plt.xlim(1, 10)
plt.ylim(-10, 5)
plt.legend()
plt.title('ERK4')
plt.show()
