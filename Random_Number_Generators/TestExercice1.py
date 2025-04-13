import unittest
from turtle import pd

from exercice1 import compute_period, compute_auto_correlation, correlation_plot, compute_period_improved_generator


class MyTestCase(unittest.TestCase):

    # Completing to generate the tables
    def test_exercice_1_period_p1(self):
        k = 899
        l = 0
        m = 32768
        x_0 = 12

        # Print the period for the initial parameters
        print('Period: ' + str(compute_period(k, l, m, x_0)))

        # Testing the variation of k
        print('k variation test:')
        for k_var in range(895, 905):
            print(k_var, compute_period(k_var, l, m, x_0))

        # Testing the variation of l
        print('l variation test:')
        for l_var in range(0, 10):
            print(l_var, compute_period(k, l_var, m, x_0))

        # Testing the variation of m
        print('m variation test:')
        for m_var in range(32764, 32774):
            print(m_var, compute_period(k, l, m_var, x_0))

        # Testing the variation of x_0
        print('x_0 variation test:')
        for x_0_var in range(8, 18):
            print(x_0_var, compute_period(k, l, m, x_0_var))


    def test_exercice_1_period_p2(self):
        k = 16807
        l = 0
        m = 6075
        x_0 = 12
        results = []

        # List of values for k
        k_values = range(16802, 16812)
        print('k variation test:')
        for k_var in range(16802, 16812):
            print(k_var, compute_period(k_var, l, m, x_0))

        # Generate results for each value of k
        results = [(k_var, compute_period(k_var, l, m, x_0)) for k_var in k_values]

        # Create LaTeX table structure
        latex_table = "\\begin{tabular}{|c|c|} \\hline\n"
        latex_table += "$k$ & $\\text{Period}$ \\\\ \\hline\n"  # Table header

        # Add rows to the table with the values of k and their corresponding periods
        for k, period in results:
            latex_table += f"${k}$ & ${period}$ \\\\ \\hline\n"

        latex_table += "\\end{tabular}"

        # Save LaTeX table to file
        with open("tableau_k_periode.tex", "w") as file:
            file.write(latex_table)

        print("LaTeX table generated in 'tableau_k_periode.tex'")


    def test_exercice_1_period_p3(self):
        k = pow(7, 5)
        l = 0
        m = pow(2, 31) - 1
        x_0 = 12
        print('Period: ' + str(compute_period(k, l, m, x_0)))

    def test_exercice_2_correlation_plot_p1(self):
        k = 899
        l = 0
        m = 32768
        x_0 = 12
        correlation_plot(k, l, m, x_0)

        # Test variations for k
        print('k variation test:')
        for k_var in range(895, 905):
            correlation_plot(k_var, l, m, x_0)

        # Test variations for l
        print('l variation test:')
        for l_var in range(0, 10):
            correlation_plot(k, l_var, m, x_0)

        # Test variations for m
        print('m variation test:')
        for m_var in range(32764, 32774):
            correlation_plot(k, l, m_var, x_0)

        # Test variations for x_0
        print('x_0 variation test:')
        for x_0_var in range(8, 18):
            correlation_plot(k, l, m, x_0_var)


    def test_exercice_2_correlation_plot_p2(self):
        k = 16807
        l = 0
        m = 6075
        x_0 = 12
        correlation_plot(k, l, m, x_0)

        # Test variations for k
        print('k variation test:')
        for k_var in range(16802, 16812):
            correlation_plot(k_var, l, m, x_0)

        # Test variations for l
        print('l variation test:')
        for l_var in range(0, 10):
            correlation_plot(k, l_var, m, x_0)

        # Test variations for m
        print('m variation test:')
        for m_var in range(6070, 6080):
            correlation_plot(k, l, m_var, x_0)

        # Test variations for x_0
        print('x_0 variation test:')
        for x_0_var in range(8, 18):
            correlation_plot(k, l, m, x_0_var)


    def test_exercice_2_correlation_plot_p3(self):
        k = pow(7, 5)
        l = 0
        m = pow(2, 31) - 1
        x_0 = 12
        correlation_plot(k, l, m, x_0)


    def test_exerice_3_autocorrelation_p1(self):
        k = 899
        l = 0
        m = 32768
        x_0 = 12

        compute_period(k, l, m, x_0)
        N = [10000, 100000]
        tau = [10, 100, 1000]

        # Compute the correlations
        correlations = []
        for n in N:
            row = []
            for t in tau:
                row.append(compute_auto_correlation(k, l, m, x_0, n, t))
            correlations.append(row)

        # Generate LaTeX table
        latex_table = "\\begin{tabular}{|c|" + "c|" * len(tau) + "} \\hline\n"
        latex_table += "$n \\backslash \\tau$ & " + " & ".join([f"${t}$" for t in tau]) + " \\\\ \\hline\n"

        for i, n in enumerate(N):
            latex_table += f"${n}$ & " + " & ".join([f"{corr:.5e}" for corr in correlations[i]]) + " \\\\ \\hline\n"

        latex_table += "\\end{tabular}"

        # Save LaTeX table to file
        with open("tableau_double_entree.tex", "w") as file:
            file.write(latex_table)

        print("LaTeX table generated in 'tableau_double_entree.tex'")


    def test_exerice_3_autocorrelation_p2(self):
        k = 16807
        l = 0
        m = 6075
        x_0 = 12

        compute_period(k, l, m, x_0)
        N = [10000, 100000]
        tau = [10, 100, 1000]

        # Compute the correlations
        correlations = []
        for n in N:
            row = []
            for t in tau:
                row.append(compute_auto_correlation(k, l, m, x_0, n, t))
            correlations.append(row)

        # Generate LaTeX table
        latex_table = "\\begin{tabular}{|c|" + "c|" * len(tau) + "} \\hline\n"
        latex_table += "$n \\backslash \\tau$ & " + " & ".join([f"${t}$" for t in tau]) + " \\\\ \\hline\n"

        for i, n in enumerate(N):
            latex_table += f"${n}$ & " + " & ".join([f"{corr:.5e}" for corr in correlations[i]]) + " \\\\ \\hline\n"

        latex_table += "\\end{tabular}"

        # Save LaTeX table to file
        with open("tableau_double_entree.tex", "w") as file:
            file.write(latex_table)

        print("LaTeX table generated in 'tableau_double_entree.tex'")


    def test_exerice_3_autocorrelation_p3(self):
        k = pow(7, 5)
        l = 0
        m = pow(2, 31) - 1
        x_0 = 12

        N = [10000, 100000]
        tau = [10, 100, 1000]

        # Compute the correlations
        correlations = []
        for n in N:
            row = []
            for t in tau:
                row.append(compute_auto_correlation(k, l, m, x_0, n, t))
            correlations.append(row)

        # Generate LaTeX table
        latex_table = "\\begin{tabular}{|c|" + "c|" * len(tau) + "} \\hline\n"
        latex_table += "$n \\backslash \\tau$ & " + " & ".join([f"${t}$" for t in tau]) + " \\\\ \\hline\n"

        for i, n in enumerate(N):
            latex_table += f"${n}$ & " + " & ".join([f"{corr:.5e}" for corr in correlations[i]]) + " \\\\ \\hline\n"

        latex_table += "\\end{tabular}"

        # Save LaTeX table to file
        with open("tableau_double_entree.tex", "w") as file:
            file.write(latex_table)

        print("LaTeX table generated in 'tableau_double_entree.tex'")


    def test_exerice_4_period_improved_generator(self):

        # k = 16807
        # l = 0
        # m = 6075
        # x_0 = 12

        k = 16807
        l = 0
        m = 6075
        x_0 = 12

        # k = 899
        # l = 1
        # m = 32768
        # x_0 = 12

        # k = pow(7, 5)
        # l = 0
        # m = pow(2, 31) - 1
        # x_0 = 12
        tau = 100
        print(compute_period_improved_generator(k, l, m, x_0, tau))
