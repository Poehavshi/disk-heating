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


if __name__ == "__main__":
    model = SumModel()
    print(model.mu_array)
