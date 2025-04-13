import matplotlib.pyplot as plt
import numpy as np

def compute_period(k, l, m, x_0):
    x = {x_0}  # create x array initialized with x0
    period_result = 1  # initialize period

    x_n = x_0  # initialize x_n
    x_n_plus_1 = (k * x_n + l) % m

    while not (x_n_plus_1 in x):  # as long as xn+1 is different from the x in the array, continue
        x.add(x_n_plus_1)
        x_n = x_n_plus_1
        x_n_plus_1 = (k * x_n + l) % m
        period_result += 1
        if period_result % 10000 == 0:
            print('current period : ' + str(period_result))

    return period_result

def correlation_plot(k, l, m, x_0):
    period_max_for_print = 100000
    x_plot = []
    y_plot = []

    period_result = 1  # initialize period

    x_n = x_0  # initialize x_n
    x_n_plus_1 = (k * x_n + l) % m

    x_plot.append(x_n)
    y_plot.append(x_n_plus_1)
    while period_result < period_max_for_print:
        x_n = x_n_plus_1
        x_n_plus_1 = (k * x_n + l) % m
        x_plot.append(x_n)
        y_plot.append(x_n_plus_1)
        period_result += 1

    print('period : ' + str(period_result))
    print(len(x_plot))
    plt.scatter(x_plot, y_plot, label='correlation', s=1)
    plt.xlabel('xn')
    plt.ylabel('xn+1')
    plt.legend()
    plt.grid(True)
    plt.show()


def compute_auto_correlation(k, l, m, x_0, N, tau):
    x = [x_0]  # create x array initialized with x0
    x_n = x_0  # initialize x_n
    x_n_plus_1 = (k * x_n + l) % m

    for i in range(1, N + tau, 1):
        x.append(x_n_plus_1)
        x_n = x_n_plus_1
        x_n_plus_1 = (k * x_n + l) % m

    correlation_numerator = 0
    correlation_denominator = 0
    x_mean = np.mean(x)

    for i in range(1, N, 1):
        correlation_numerator += (x[i] - x_mean) * (x[i + tau] - x_mean)
        correlation_denominator += pow((x[i] - x_mean), 2)

    correlation = correlation_numerator / correlation_denominator
    return correlation


def compute_period_improved_generator(k, l, m, x_0, tau):
    k_1 = 899
    l_1 = 0
    m_1 = 32768
    x_0_1 = 12

    m1_size = 800
    m1 = []
    N = 100000  # number of elements for computing correlation
    N_mean = N - tau - 1

    x_n_1 = x_0_1
    x_n_plus_1_1 = (k_1 * x_n_1 + l_1) % m_1

    for i in range(0, m1_size):
        m1.append(x_n_plus_1_1)
        x_n_1 = x_n_plus_1_1
        x_n_plus_1_1 = (k_1 * x_n_1 + l_1) % m_1

    # calculate the index with the second generator
    x_n = x_0
    x_n_plus_1 = (k * x_n + l) % m
    index = int(x_n_plus_1 * m1_size / m)
    x_n = x_n_plus_1

    x_n_plus_1_improved = m1[index]

    # replace with a number drawn from the first generator
    x_n_plus_1_1 = (k_1 * x_n_1 + l_1) % m_1
    m1[index] = x_n_plus_1_1
    x_n_1 = x_n_plus_1_1

    x = []  # create x array initialized with x0
    period_result = 1  # initialize period
    # improved generator

    for i in range(1, N, 1):  # as long as xn+1 is different from the x in the array, continue
        x.append(x_n_plus_1_improved)

        x_n_plus_1 = (k * x_n + l) % m
        index = int(x_n_plus_1 * m1_size / m)
        x_n = x_n_plus_1

        x_n_plus_1_improved = m1[index]

        x_n_plus_1_1 = (k_1 * x_n_1 + l_1) % m_1
        m1[index] = x_n_plus_1_1
        x_n_1 = x_n_plus_1_1

        period_result += 1
        # print(x_n_plus_1_improved)

    correlation_numerator = 0
    correlation_denominator = 0
    x_mean = np.mean(x)

    for i in range(1, N_mean, 1):
        correlation_numerator += (x[i] - x_mean) * (x[i + tau] - x_mean)
        correlation_denominator += pow((x[i] - x_mean), 2)

    correlation = correlation_numerator / correlation_denominator
    return correlation
