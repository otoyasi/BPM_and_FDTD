
# -*- coding: utf-8 -*-
"""
交互方向陰解法による2次元の熱拡散方程式の計算
"""

import matplotlib.animation as animation  # アニメーション作成のためのメソッドをインポート
import pml
import par
import neff
import index
from sklearn.metrics import completeness_score
import numpy as np
import numpy.matlib
from decimal import Decimal, ROUND_HALF_UP
import scipy.sparse.linalg as spla
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import time
import math
np.set_printoptions(threshold=np.inf)


dtnum = 301  # 時間をいくつに分割して計算するか 一の位を1とした方が良い 時刻tのfor文の終わりの数字に関係あり
dxnum = 175  # xとyをいくつに分割して計算するか
dynum = 175

thick = 10  # x方向の大きさ
width = 10  # y方向の大きさ
time_calc = 500  # 計算する時間
beta = 0.1  # 上の式中ではラムダ　温度拡散率 lambda = 熱伝導率/(密度*比熱)

"""計算準備"""
# 空の解を用意
solution_k1 = np.zeros([dxnum, dynum], dtype=np.complex128)
# solution_k2の初期値が初期条件になる
solution_k2 = np.zeros([dxnum, dynum], dtype=np.complex128)


# 境界条件
irr_boundary = 1e-30  # 表面の境界条件 (0,t)における温度
rear_boundary = 1e-30  # 裏面と反対側の境界条件　(x,t)における温度 0を基準となる温度とした
side_0 = 1e-30  # y=0の温度
side_y = 1e-30  # y=0とは反対側の境界の温度

dt = time_calc / (dtnum - 1)
dx = thick / (dxnum - 1)
dy = width / (dynum - 1)
B = (2*((dx/dy)**2)-((dx**2)/(beta*dt)))

e = -(2+(dy**2)/(beta*dt))
C = (2*((dy/dx)**2)-((dy**2)/(beta*dt)))


d = np.zeros([10, 10])
a1 = np.zeros([dxnum, dynum], dtype=np.complex128)
c1 = np.zeros([dxnum, dynum], dtype=np.complex128)
b1 = np.zeros([dxnum, dynum], dtype=np.complex128)
d1 = np.zeros([dxnum, dynum], dtype=np.complex128)

a2 = np.zeros([dxnum, dynum], dtype=np.complex128)
c2 = np.zeros([dxnum, dynum], dtype=np.complex128)
b2 = np.zeros([dxnum, dynum], dtype=np.complex128)
d2 = np.zeros([dxnum, dynum], dtype=np.complex128)


file_index = open("index.csv", "w")
file_matrix= open("matrix.csv", "w")
file_arry = open("arry.csv","w")
be = neff.value*par.w0
IM = 4j*be/par.dz
print("wo=",par.w0)
for i in range(dxnum):
    for j in range(dynum):
        #print("j= ",j)
        a1[i, j] = -(1 / (par.dx**2))
        c1[i, j] = -(1 / (par.dx**2))
        b1[i, j] = (4j*par.w0*(index.n[i, j]**2))/((par.c**2)*par.dt)  \
                     + 2/(par.dz**2)-((par.w0**2) * (index.n[i, j]**2))/(2*(par.c**2))
        #b1[i,j]=  (index.n[i, j]**2)/(2*(par.c**2))
        #b1[i,j]= 1e-10
        #print("b1[i,j]=",b1[i,j])
    for j in range(dynum):
        # print("d[i,j=",d[i,j])
        a_matrix_1 = np.identity(dxnum) *b1[i, j] \
            + np.eye(dxnum, k=-1)*a1[i, j] \
            + np.eye(dxnum, k=1)*c1[i, j]
        # print("a_matrix=\n",a_matrix_1)


"""Ax=bのAを用意"""
for i in range(dxnum):

    for j in range(dynum):
        #print("j= ",j)
        a2[i, j] = -(1 / (par.dx**2))
        c2[i, j] = -(1 / (par.dx**2))
        b2[i, j] = (4j*par.w0*(index.n[i, j]**2))/((par.c**2)*par.dt)  \
                    + 2/(par.dz**2)-((par.w0**2) * (index.n[i, j]**2))/(2*(par.c**2))

    for j in range(dynum):
        file_index.write("   "+str(j)+",")
        file_index.write("   "+str(i)+",")
        file_index.write("   "+str((b2[i,j]))+"\n")
#file.write("   "+str(solution.imag[:,:,k])+"\n")

        a_matrix_2 = np.identity(dynum) * b2[i, j] \
            + np.eye(dynum, k=-1)*a2[i, j]  \
            + np.eye(dynum, k=1)*c2[i, j]
        # print("a_matrix_2=\n",a_matrix_2)

        if j == dynum:
            file_index.write("\n")

file_index.write("   "+str(j)+",")
file_index.write("   "+str(i)+",")
file_index.close()


# 疎行列を格納　csr方式
a_matrix_csr1 = csr_matrix(a_matrix_1)
a_matrix_csr2 = csr_matrix(a_matrix_2)
file_matrix.write("  "+str(a_matrix_1)+"\n")
file_matrix.write("  "+str(a_matrix_2)+"\n")
file_matrix.close()
# print("barray=",a_matrix_csr1)

# ADI法ではk+1時刻とk+2時刻を計算するので、for文の回数は半分になる
number = Decimal(str(dtnum/2)).quantize(Decimal('0'), rounding=ROUND_HALF_UP)
number = int(number)

solution = np.zeros([dxnum, dynum, number], dtype=np.complex128)  # 解を代入する行列を作成

# bpmは、ずれた列を作成して、足すためのもの
bpm1 = np.zeros([dxnum, dynum+2], dtype=np.complex128)
bpm1[:, 0] = side_0
bpm1[:, -1] = side_y

# for k in range(number):
#     for i in range(dxnum):
#         for j in range(dynum):
solution_k2[88, 88] = 1e1
            #bpm1[i, j] = 1

bpm2 = np.zeros([dxnum+2, dynum], dtype=np.complex128)
bpm2[0, :] = irr_boundary
bpm2[-1, :] = rear_boundary

# ADI法の計算
for k in range(number):  # 時刻tについて
    print("k=", k)
    for j in range(dynum-1):  # xを計算し、yの数繰り返すことで、時刻k+1のTijを計算
        bpm1[:, 1:dynum+1] = solution_k2
        b_array_1 = b1[:,j+1] * bpm1[:, j+1] \
            + a1[:,j] * (bpm1[:, j] + bpm1[:, j+2])

        # print("k=",k)
        # print("j=",j)
        #print("B=",b_array_1)
        #print("B=",b1[k][j])
        # bの最初と最後に境界条件
        # b_array_1[0] -= irr_boundary
        # b_array_1[-1] -= rear_boundary

        # 解を求める
        bpm_solve1 = spla.dsolve.spsolve(a_matrix_csr1, b_array_1)  # xについての解
        solution_k1[:, j] = bpm_solve1

    for i in range(dxnum-1):  # yを計算し、xの数繰り返すことで、時刻k+2のTijを計算
        bpm2[1:dxnum+1, :] = solution_k1
        b_array_2 = b2[i+1,:] * bpm2[i+1, :] \
            + a2[i,:] * (bpm2[i, :]) -c2[i,:] *(bpm2[i+2, :])
        file_arry.write("   "+str(b_array_2)+"\n")
        # bの最初と最後に境界条件
        # b_array_2[0] -= side_0
        # b_array_2[-1] -= side_y

        # 解を求める
        bpm_solve2 = spla.dsolve.spsolve(a_matrix_csr2, b_array_2)  # yについての解
        solution_k2[i, :] = bpm_solve2
    solution[:, :, k] = solution_k2
file_arry.close()

file = open("solution_imag.csv", "w")
for i in range(dynum):  # xを計算し、yの数繰り返すことで、時刻k+1のTijを計算
    if i % 4 == 0:
        for j in range(dxnum):  # yを計算し、xの数繰り返すことで、時刻k+2のTijを計算

            file.write("   "+str(i)+",")
            file.write("   "+str(j)+",")
            file.write("   "+str(math.sqrt(solution.real[j,i,90]**2+solution.imag[j, i,90]**2))+"\n")
# #file.write("   "+str(solution.imag[:,:,k])+"\n")
#         if i == 153:
#            file.write("\n")
# for k in range(1,number):
#     for j in range(1, dynum):  # xを計算し、yの数繰り返すことで、時刻k+1のTijを計算
#         for i in range(1, dxnum):  # yを計算し、xの数繰り返すことで、時刻k+2のTijを計算

#             file.write("   "+str(k)+",")
#             file.write("   "+str(j)+",")
#             file.write("   "+str(i)+",")
#             file.write("   "+str((solution.imag[j,i,k])**2)+"\n")

# #file.write("   "+str(solution.imag[:,:,k])+"\n")
#         if i == 153:
#            file.write("\n")
#file.write("   "+str(solution.imag[:,25:150,11])+"\n")
# file.write(str(solution.imag[:,:,100]**2))


# for k in range(1, number):
#     if k % 5 == 0:
#         for i in range(25, 154):  # yを計算し、xの数繰り返すことで、時刻k+2のTijを計算

#             file.write("   "+str(k)+",")
#             #file.write("   "+str(j)+",")
#             file.write("   "+str(i)+",")
#             file.write("   "+str((solution.imag[10, i, k])**2)+"\n")
# #file.write("   "+str(solution.imag[:,:,k])+"\n")
#             if i == 153:
#                 file.write("\n")
file.close()

ax = sns.heatmap(solution.imag[:, :, 10], linewidth=0, vmin=0, vmax=100)
plt.show()

ax = sns.heatmap(solution.imag[:, :, 50], linewidth=0, vmin=0, vmax=100)
plt.show()

# ax = sns.heatmap(solution.imag[:, :, 100], linewidth=0, vmin=0, vmax=100)
# plt.show()

# ax = sns.heatmap(solution.imag[:, :, 1000], linewidth=0, vmin=0, vmax=100)
# plt.show()

# ax = sns.heatmap(solution.imag[:, :, 2000], linewidth=0, vmin=0, vmax=100)
# plt.show()

# ax = sns.heatmap(solution.imag[:, 10, :], linewidth=0, vmin=0, vmax=100)
# plt.show()

# ax = sns.heatmap(solution.imag[:, :, 500], linewidth=0, vmin=0, vmax=100)
# plt.show()

# ax = sns.heatmap(solution.imag[:, :, 2000], linewidth=0, vmin=0, vmax=100)
# plt.show()


fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(projection='3d')

x = list(range(dxnum))
y = list(range(dynum))
X, Y = np.meshgrid(x, y)


def update(i, solution, fig_title):
    if i != 0:
        ax.clear()
    ax.view_init(elev=60., azim=60.)  # アングル設定
    # ax.plot_surface(X,Y,u[X,Y,i],rstride=1,cstride=1,cmap='Greys',linewidth=0.3) # サーフェスプロット
    ax.plot_wireframe(X, Y, solution[X, Y, i], rstride=1,
                      cstride=1, color='blue', linewidth=0.3)  # ワイヤーフレーム表示

    ax.set_title(fig_title + 'i=' + str(i))
    ax.set_zlim()  # zの描画範囲の指定
    ax.set_xlabel('X')  # ラベル
    ax.set_ylabel('Y')
    ax.set_zlabel('U')

    return ax


ani = animation.FuncAnimation(fig, update, fargs=(
    solution.imag, 'Wave motion: time step='), interval=1, frames=number, blit=False)
fig.show()

ani.save("output.gif", writer="pillow")
