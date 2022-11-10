import tkinter as tk

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
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


class Program:
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
        self.root.resizable(width=False, height=False)
        self.fig = plt.figure()

        self.canvas = FigureCanvasTkAgg(self.fig, master=self.root)
        self.canvas.get_tk_widget().grid(rowspan=20, column=1)

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

        self.entry_I_value = tk.Entry(textvariable=tk.StringVar(self.root, value='1000'))
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

        self.current_variant = 0
        self.var = tk.IntVar()
        self.var.set(0)
        self.r_variant = tk.Radiobutton(text="r", command=self.select_another_variable, variable=self.var, value=0)
        self.r_variant.grid(row=7, column=2, padx=30)

        self.t_variant = tk.Radiobutton(text="t", command=self.select_another_variable, variable=self.var, value=1)
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

        self.diff_button = tk.Button(text='find diff', command=self.find_diff)
        self.diff_button.grid(row=11, column=2, columnspan=2)

        self.root.mainloop()

    def select_another_variable(self):
        self.current_variant = self.var.get()

    def generate_matrix(self, hr_coef=1.0, ht_coef=1.0):
        # шаги сетки
        self.I = self.entry_I_value.get()
        self.R = self.entry_R_value.get()
        self.hr = (float(self.R) / float(self.I)) * hr_coef

        self.K = self.entry_K_value.get()
        self.T = self.entry_T_value.get()
        self.ht = (float(self.T) / float(self.K)) * ht_coef

        self.mas_r = [self.hr * i for i in range(int(self.I) + 1)]
        self.mas_t = [self.ht * i for i in range(int(self.K) + 1)]
        self.matrix = [[0 for i in range(int(self.I) + 1)] for j in range(int(self.K) + 1)]

        for i in range(len(self.matrix[0])):
            self.matrix[0][i] = self.phi_r(self.mas_r[i])
        for i in range(1, len(self.matrix)):
            self.matrix[i][-1] = self.Ub

        for i in range(1, len(self.matrix)):
            self.E = 1 + 2 * self.alpha * self.ht / self.c / self.l + 4 * self.k * self.ht / self.c / ((self.hr) ** 2)
            self.F = - 4 * self.k * self.ht / self.c / ((self.hr) ** 2)
            self.ps = [- self.F / self.E]
            self.qs = [self.matrix[i - 1][0] / self.E]
            # self.qs=[self.matrix[i-1][0]/(1 + 2 * ( self.alpha * self.ht / ( self.c * self.l ) ) + 4 * ( self.k * self.ht / ( self.c * self.hr * self.hr ) ))]
            # self.ps=[-(-4 * self.k * self.ht / ( self.c * self.hr * self.hr )) / (1 + 2 * ( self.alpha * self.ht / ( self.c * self.l ) ) + 4 * ( self.k * self.ht / ( self.c * self.hr * self.hr ) ))]
            for j in range(1, len(self.matrix[i]) - 1):
                self.A = 1 + 2 * self.alpha * self.ht / self.c / self.l + 2 * self.k * self.ht / self.c / (
                        (self.hr) ** 2)
                self.B = - (self.k * self.ht / self.hr / self.c) * (
                        (2 * self.mas_r[j] + self.hr) / (2 * self.mas_r[j] * self.hr))
                self.C = - (self.k * self.ht / self.hr / self.c) * (
                        (2 * self.mas_r[j] - self.hr) / (2 * self.mas_r[j] * self.hr))
                self.ps.append(-self.B / (self.C * self.ps[j - 1] + self.A))
                self.qs.append((self.matrix[i - 1][j] - self.C * self.qs[j - 1]) / (self.C * self.ps[j - 1] + self.A))
                # self.qs.append(( self.matrix[i-1][j] - self.qs[j-1] * ( -self.k * self.ht * ( 1/(self.hr*self.hr) - 1/(2*self.hr*self.mas_r[j]) ) / self.c ) ) / ( ( -self.k * self.ht * ( 1/(self.hr*self.hr) - 1/(2*self.hr*self.mas_r[j]) ) / self.c ) * self.ps[j-1] + (1/self.ht + 2*(self.alpha/(self.c*self.l)) + 2*self.k/(self.c*self.hr*self.hr) ) ))
                # self.ps.append(-(-self.k*self.ht*(1/(2*self.hr*self.mas_r[j])+1/(self.hr*self.hr))/self.c) / ( (( -self.k * self.ht * ( 1/(self.hr*self.hr) - 1/(2*self.hr*self.mas_r[j]) ) / self.c ) * self.ps[j-1]) +  (1/self.ht + 2*(self.alpha/(self.c*self.l)) + 2*self.k/(self.c*self.hr*self.hr) )) )
            for j in range(len(self.matrix[i]) - 2, -1, -1):
                if j == len(self.matrix[i]) - 2:
                    self.matrix[i][j] = self.qs[j]
                elif j != len(self.matrix[i]) - 2:
                    self.matrix[i][j] = self.matrix[i][j + 1] * self.ps[j] + self.qs[j]
                #    self.matrix[i][j] = self.matrix[i-1][j] + self.ht * (-2 * ( self.alpha / (self.c * self.l) ) * (self.matrix[i-1][j] - self.Uc ) - 4 * ( (self.k * ( self.matrix[i-1][j] - self.matrix[i-1][j+1] ) ) / (self.c * self.hr * self.hr ) ) )
                ##elif j==1:
                ##    self.matrix[i][j] = ( (( self.matrix[i][j-1] - self.matrix[i-1][j-1] ) / self.ht ) + ( ( 2 * self.alpha * self.matrix[i][j-1] ) / ( self.c * self.l ) ) ) * ( ( self.hr * self.hr ) / ( 2 * ( self.k / ( self.c ) ) ) ) + 2 * self.matrix[i][j-1]
                ##else:
                ##    self.matrix[i][j] = ( (( self.matrix[i][j-1] - self.matrix[i-1][j-1] ) / self.ht ) + ( ( 2 * self.alpha * self.matrix[i][j-1] ) / ( self.c * self.l ) ) - ( self.k / ( self.c * self.mas_r[j] ) ) * ( - ( ( self.matrix[i][j-2] ) / ( 2 * self.hr ) ) + self.mas_r[j] * ( ( -2 * self.matrix[i][j-1] + self.matrix[i][j-2] ) / ( self.hr * self.hr ) ) ) ) / ( ( self.k * ( self.hr - 2 * self.mas_r[j] ) ) / ( self.c * self.mas_r[j] * 2 * self.hr ) )
                # else:
                #    pass

    def build(self):
        self.generate_matrix()
        if self.current_variant == 0:  # если выбрано r
            ind = 0
            while self.mas_r[ind] < float(self.entry_value.get()):
                ind += 1
            plt.plot(self.mas_t, [self.matrix[i][ind] for i in range(len(self.matrix))],
                     label=f"I={self.I}, K={self.K}")
        else:
            ind = 0
            while self.mas_t[ind] < float(self.entry_value.get()):
                ind += 1
            plt.plot(self.mas_r, [self.matrix[ind][i] for i in range(len(self.matrix[ind]))],
                     label=f"I={self.I}, K={self.K}")
        plt.legend()
        plt.gcf().canvas.draw()
        self.canvas.draw()

    def build_analytical(self):
        R = self.entry_R_value.get()
        T = self.entry_T_value.get()

        if self.current_variant == 0:  # если выбрано r
            summodel = SumModel()
            plot_data = summodel.generate_w_data(100, float(self.entry_value.get()), "r", int(T))
            plt.plot(plot_data[0], plot_data[1], label="analytical")
        if self.current_variant == 1:  # если выбрано t
            summodel = SumModel()
            plot_data = summodel.generate_w_data(100, float(self.entry_value.get()), "t", int(R))
            plt.plot(plot_data[0], plot_data[1], label="analytical")
        plt.gcf().canvas.draw()
        self.canvas.draw()

    def clear(self):
        plt.clf()
        self.canvas.draw()

    def find_diff(self):
        self.generate_matrix()
        summodel = SumModel()

        max_eps = 0
        max_pair = None
        # for i in range(int(self.I)):
        #     for k in range(1, int(self.K)):
        i = int(self.I) // 2
        k = int(self.K) // 2

        precise_solve = summodel.calculate_sum(self.mas_r[i], self.mas_t[k], 100)
        unprecise_solve = self.matrix[k][i]
        eps = abs(precise_solve - unprecise_solve)
        if eps > max_eps:
            max_eps = eps
            max_pair = (k,i)
        print("eps", max_eps)
        print("max pair", max_pair)


    def phi_r(self, value):
        if float(value) <= float(self.entry_R_value.get()) / 4:
            return 10
        return 0


if __name__ == '__main__':
    main = Program()
