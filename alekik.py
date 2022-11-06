import tkinter as tk

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

"""
Модель - представление данных для программы
Содержит:
* generate_exp_data - функция, которая генерирует данные вида (x, y), подходящие для отображения в matplotlib.
"""

# !TODO расчёт суммы ряда

import numpy as np
from numpy import exp
from scipy.special import jv, jn_zeros


def generate_exp_data(base, exponent):
    x = np.linspace(0, 10, 800)
    y = (x * base) ** exponent
    return x, y


class SumModel:
    def __init__(self,
                 R=4,
                 l=0.5,
                 u_c=0,
                 u_b=0,
                 alpha=0.005,
                 T=150,
                 k=0.065,
                 c=1.35
                 ):
        self.c = c
        self.k = k
        self.T = T
        self.alpha = alpha
        self.u_b = u_b
        self.u_c = u_c
        self.l = l
        self.R = R

        self.mu_array = jn_zeros(0, 1251)

    def _calculate_term(self, n: int, r: float, t: float) -> float:
        """
        Функция подсчёта n-ого слагаемого суммы

        :param n: порядок слагаемого.
        :param r: аргумент функции w.
        :param t: аргумент функции w.

        :return значение одного слагаемого суммы
        """
        # mu_n = jn_zeros(0, n)[n - 1]
        mu_n = self.mu_array[n - 1]
        result = (5 * jv(1, mu_n / 4)) / (mu_n * (jv(1, mu_n)) ** 2)
        # result=(5 * jv(1, mu_n / 4)) / (mu_n )
        result *= exp(-(t * (self.l * self.k * (mu_n / self.R) ** 2 + 2 * self.alpha)) / (self.l * self.c))
        result *= jv(0, (mu_n * r) / self.R)
        return result

    def phi(self, N, t):
        result = ((self.R ** 2) * self.c * 5 * 2 ** 0.5) / (2 * self.k * np.pi ** 2 * t * (N - 0.25) ** 1.5)
        result *= np.exp((-(2 * self.alpha * t) / (self.l * self.c)) - (
                (t * self.k * (np.pi ** 2) * ((N - 0.25) ** 2)) / (self.c * (self.R ** 2))))

        return result

    def calculate_number_of_iterations(self, epsilon, t):
        N = 1
        while self.phi(N, t) > epsilon:
            N += 1
        return N

    def calculate_eps_of_iterations(self, N, t):
        return self.phi(N, t)

    def calculate_sum(self, r: float, t: float, N: int) -> float:
        """
        Функция подсчёта значения функции w(r, t)

        :param r: аргумент функции w.
        :param t: аргумент функции w.
        :param N: количество элементов ряда

        :return значение функции w(r,t)
        """
        result = 0
        for i in range(N):
            result += self._calculate_term(i + 1, r, t)
        return result

    def generate_w_data(self, N: int, r: float, p: str, x: int):
        """
        Генерирует значения функции w(r, t)

        :param N: количество элементов (точность подсчёта функции)
        :param r: фиксированный параметр.
        :param p: строка с тем параметром, который является фиксированным.
        :param x: максимальное значение изменяемого параметра

        :return вектор двух numpy массивов
        """
        ox = np.linspace(0.001, x, 800)
        w = np.zeros(800)
        if p == 'r':
            for i in range(800):
                w[i] = self.calculate_sum(r=r, t=ox[i], N=N)
        else:

            for i in range(800):
                w[i] = self.calculate_sum(r=ox[i], t=r, N=N)
        if r != 0:
            # print(abs(abs(self.calculate_sum(r=r,t=ox[0],N=N))-abs(self.calculate_sum(r=r,t=ox[0],N=self.calculate_number_of_iterations(E,r)))))
            print(w[0])
        return ox, w


class programm():
    def __init__(self):

        self.alpha = 0.005
        self.l = 0.5
        self.Uc = 0
        self.Ub = 0
        self.c = 1.35
        self.k = 0.065

        self.matrix = []
        self.ht = 0
        self.hr = 0
        self.root = tk.Tk()
        self.root.geometry("960x480")
        self.root.resizable(width=0, height=0)
        self.fig = plt.figure()

        self.canv = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canv.get_tk_widget().grid(rowspan=20, column=1)

        # R
        self.label_R_value = tk.Label(text="Input the R value")
        self.label_R_value.grid(row=2, column=2, padx=30)

        self.entry_R_value = tk.Entry(textvariable=tk.StringVar(self.root, value='4'))
        self.entry_R_value.grid(row=2, column=3, padx=40)

        # T
        self.label_T_value = tk.Label(text="Input the T value")
        self.label_T_value.grid(row=3, column=2, padx=30)

        self.entry_T_value = tk.Entry(textvariable=tk.StringVar(self.root, value='150'))
        self.entry_T_value.grid(row=3, column=3, padx=40)

        # I
        self.label_I_value = tk.Label(text="Input the I value")
        self.label_I_value.grid(row=4, column=2, padx=30)

        self.entry_I_value = tk.Entry(textvariable=tk.StringVar(self.root, value='20'))
        self.entry_I_value.grid(row=4, column=3, padx=40)

        # K
        self.label_K_value = tk.Label(text="Input the K value")
        self.label_K_value.grid(row=5, column=2, padx=30)

        self.entry_K_value = tk.Entry(textvariable=tk.StringVar(self.root, value='1000'))
        self.entry_K_value.grid(row=5, column=3, padx=40)

        # value
        self.label_value = tk.Label(text="Input the select value")
        self.label_value.grid(row=1, column=2, padx=30)

        self.entry_value = tk.Entry(textvariable=tk.StringVar(self.root, value='0'))
        self.entry_value.grid(row=1, column=3, padx=40)

        # �����
        self.current_variant = 0
        self.var = tk.IntVar()
        self.var.set(0)
        self.r_variant = tk.Radiobutton(text="r", command=self.change, variable=self.var, value=0)
        self.r_variant.grid(row=7, column=2, padx=30)

        self.t_variant = tk.Radiobutton(text="t", command=self.change, variable=self.var, value=1)
        self.t_variant.grid(row=7, column=3, padx=30)

        # button build
        self.build_button = tk.Button(text="run", command=self.build)
        self.build_button.grid(row=8, column=2, columnspan=2)

        # button analytical build
        self.a_build_button = tk.Button(text="analytical build", command=self.build_analytical)
        self.a_build_button.grid(row=9, column=2, columnspan=2)

        # button clear
        self.clear_button = tk.Button(text="clear", command=self.clear)
        self.clear_button.grid(row=10, column=2, columnspan=2)

        self.root.mainloop()

    def change(self):
        self.current_variant = self.var.get()

    def generate_matrix(self):
        self.hr = float(self.entry_R_value.get()) / float(self.entry_I_value.get())
        self.ht = float(self.entry_T_value.get()) / float(self.entry_K_value.get())
        self.mas_x = [self.hr * i for i in range(int(self.entry_I_value.get()) + 1)]
        self.mas_y = [self.ht * i for i in range(int(self.entry_K_value.get()) + 1)]
        self.matrix = [[0 for i in range(int(self.entry_I_value.get()) + 1)] for j in
                       range(int(self.entry_K_value.get()) + 1)]

        for i in range(len(self.matrix[0])):
            self.matrix[0][i] = self.phi_r(self.mas_x[i])
        for i in range(1, len(self.matrix)):
            self.matrix[i][-1] = self.Ub

        for i in range(1, len(self.matrix)):
            for j in range(len(self.matrix[i]) - 1):
                if j == 0:
                    self.matrix[i][j] = self.matrix[i - 1][j] + self.ht * (
                            -2 * (self.alpha / (self.c * self.l)) * (self.matrix[i - 1][j] - self.Uc) - 4 * (
                            (self.k * (self.matrix[i - 1][j] - self.matrix[i - 1][j + 1])) / (
                            self.c * self.hr * self.hr)))
                else:
                    self.matrix[i][j] = self.matrix[i - 1][j] + self.ht * (
                            -2 * (self.alpha / (self.c * self.l)) * (self.matrix[i - 1][j] - self.Uc) + (
                            self.k / (self.c * self.mas_x[j])) * (
                                    ((self.matrix[i - 1][j + 1] - self.matrix[i - 1][j - 1]) / (2 * self.hr)) +
                                    self.mas_x[j] * ((self.matrix[i - 1][j + 1] - 2 * self.matrix[i - 1][j] +
                                                      self.matrix[i - 1][j - 1]) / (
                                                             self.c * self.hr * self.hr))))

    def build(self):
        self.generate_matrix()
        if (self.current_variant == 0):
            ind = 0
            while (self.mas_x[ind] < float(self.entry_value.get())):
                ind += 1
            # print(ind,self.mas_x)
            # print([self.matrix[i][ind] for i in range(len(self.matrix))])
            # print(self.hr)
            plt.plot(self.mas_y, [self.matrix[i][ind] for i in range(len(self.matrix))],
                     label=f"I={self.entry_I_value.get()}, K={self.entry_K_value.get()}")
        else:
            ind = 0
            while (self.mas_y[ind] < float(self.entry_value.get())):
                ind += 1
            # print(ind,self.mas_y)
            # print([self.matrix[ind][i] for i in range(len(self.matrix[ind]))])
            plt.plot(self.mas_x, [self.matrix[ind][i] for i in range(len(self.matrix[ind]))],
                     label=f"I={self.entry_I_value.get()}, K={self.entry_K_value.get()}")
        plt.gcf().canvas.draw()
        self.canv.draw()

    def build_analytical(self):
        R = self.entry_R_value.get()
        T = self.entry_T_value.get()

        if self.current_variant == 0:  # если выбрано r
            summodel = SumModel()
            plot_data = summodel.generate_w_data(100, float(self.entry_value.get()), "r", int(T))
            plt.plot(plot_data[0], plot_data[1], label=f"analytical r={self.entry_value.get()}")
        if self.current_variant == 1:  # если выбрано t
            summodel = SumModel()
            plot_data = summodel.generate_w_data(100, float(self.entry_value.get()), "t", int(R))
            plt.plot(plot_data[0], plot_data[1], label=f"analytical t={self.entry_value.get()}")
        plt.gcf().canvas.draw()
        self.canvas.draw()

    def clear(self):
        plt.clf()
        self.canvas.draw()

    def phi_r(self, value):
        if (float(value) <= float(self.entry_R_value.get()) / 4):
            return 10
        return 0


main = programm()
