import tkinter as tk
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


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

        self.root.mainloop()

    def select_another_variable(self):
        self.current_variant = self.var.get()

    def generate_matrix(self):
        # шаги сетки
        I = self.entry_I_value.get()
        R = self.entry_R_value.get()
        self.hr = float(R) / float(I)

        K = self.entry_K_value.get()
        T = self.entry_T_value.get()
        self.ht = float(T) / float(K)

        # Создаём массив значений для r и t
        self.mas_r = [self.hr * i for i in range(int(I) + 1)]
        self.mas_t = [self.ht * i for i in range(int(K) + 1)]

        # список из нулей
        self.matrix = [[0 for _ in range(int(I) + 1)] for _ in range(int(K) + 1)]

        # граничные условия
        # v_i ^ 0 = ψ_i, i =¯(0, I);
        for i in range(len(self.matrix[0])):
            self.matrix[0][i] = self.phi_r(self.mas_r[i])

        # v_I^k=u_b,k=¯(1,K);
        for k in range(1, len(self.matrix)):
            self.matrix[k][-1] = self.Ub

        # todo написать прямой ход вычисления прогоночных коэффициентов
        for k in range(1, len(self.matrix)):
            for i in range(len(self.matrix[k]) - 1):
                A = self.alpha / (self.c * self.l)
                C = 2 * self.k / (self.c * self.hr ** 2)

                delta = C * self.ht
                mu = 1 + A * self.ht + C * self.ht
                if i == 0:
                    p0 = delta / mu
                    q0 = delta / mu * self.matrix[k-1][1] + self.matrix[k-1][0]
                    self.matrix[k][0] = p0 * self.matrix[k][1] + q0
                else:
                    self.matrix[k][i] = self.matrix[k - 1][i] + self.ht * (
                            -2 * (self.alpha / (self.c * self.l)) * (self.matrix[k - 1][i] - self.Uc) + (
                            self.k / (self.c * self.mas_r[i])) * (
                                    ((self.matrix[k - 1][i + 1] - self.matrix[k - 1][i - 1]) / (2 * self.hr)) +
                                    self.mas_r[i] * ((self.matrix[k - 1][i + 1] - 2 * self.matrix[k - 1][i] +
                                                      self.matrix[k - 1][i - 1]) / (
                                                             self.c * self.hr * self.hr))))
        # todo а тут обратный ход

    def build(self):
        self.generate_matrix()
        if self.current_variant == 0:
            ind = 0
            while self.mas_r[ind] < float(self.entry_value.get()):
                ind += 1
            # print(ind,self.mas_x)
            # print([self.matrix[i][ind] for i in range(len(self.matrix))])
            # print(self.hr)
            plt.plot(self.mas_t, [self.matrix[i][ind] for i in range(len(self.matrix))])
        else:
            ind = 0
            while self.mas_t[ind] < float(self.entry_value.get()):
                ind += 1
            # print(ind,self.mas_y)
            # print([self.matrix[ind][i] for i in range(len(self.matrix[ind]))])
            plt.plot(self.mas_r, [self.matrix[ind][i] for i in range(len(self.matrix[ind]))])
        plt.gcf().canvas.draw()
        self.canvas.draw()

    def phi_r(self, value):
        if float(value) <= float(self.entry_R_value.get()) / 4:
            return 10
        return 0


if __name__ == '__main__':
    main = Program()
