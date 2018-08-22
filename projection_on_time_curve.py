#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 2017
@author: Tim Veenstra
"""

import sys, getopt
import sympy as sp
import numpy as np
import scipy.odr as odr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
x, y, Xi, Yi = sp.symbols("x, y, Xi, Yi")


""" ---------------------------- PARAMETERS -------------------------------"""

knockdown_input_file = ""
time_series_input_file = ""
output_filename = ""


help_text = "projection_on_time_curve.py, Arguments:\n -h\n"
help_text += " -t   <input filename time series>\n"
help_text += " -d   <input filename gene knockdown data>\n"
help_text += " -o   <filename output txt file> \n"
help_text += " -f   <fit polynomial order, integer 1..6> \n"
help_text += " [-p] <filename output pdf plot> \n\n"
help_text += "General format input files:\n  -Tab seperated .txt files\n  -Comments and headers should start with a '#'\n"
help_text += "Format time series input file: \n  #Day \t PC1 \t PC2\n"
help_text += "Format gene knockdown data input file: \n  #Gene \t PC1 \t PC2\n"
help_text += "The order of the columns is important, the header is not\n\n"

savefig = False


def main(argv):
    global time_series_input_file, knockdown_input_file, fit_polynom_order,\
            output_filename, savefig, inputfile_figure
    try:
        opts, args = getopt.getopt(argv, "ht:d:o:f:p:")
    except getopt.GetoptError:
        print(help_text)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print(help_text)
            sys.exit()
        elif opt == "-t":
            time_series_input_file = arg
        elif opt == "-d":
            knockdown_input_file = arg
        elif opt == "-o":
            output_filename = arg
        elif opt == "-f":
            fit_polynom_order = int(arg)
        elif opt == "-p":
            savefig = True
            inputfile_figure = arg


main(sys.argv[1:])

""" ----------------------------------------------------------------------"""

"""
    Reads the data from the inputfiles
"""
try:
    data = np.genfromtxt(knockdown_input_file, delimiter='\t', usecols=[0, 1, 2],
                         dtype=None, names=['gene', 'PC1', 'PC2'])
    time_series = np.genfromtxt(time_series_input_file, delimiter='\t')
except IOError:
    print("\nERROR: failed to open input files. Make sure you have entered a correct filepath.\n")
    print(help_text)
    sys.exit()


(gene, data_x, data_y) = (data['gene'], data['PC1'], data['PC2'])

(day, time_x, time_y) = (time_series[:, 0], time_series[:, 1],
                         time_series[:, 2])

""" ---------------------------------------------------------------------- """

"""
    Fits a function to the timeseries
    (For the documentation of the odr-module, see:
        https://docs.scipy.org/doc/scipy/reference/odr.html)
"""


def func(B, x):
    """
        Function for the polynomial. B is an array containing the parameters
        for the polynomial. Returns B[0] + B[1]*x + ... + B[n]*x^n where
        n is the specified order of the polynomial.
    """
    global fit_polynom_order
    return sum(([B[i]*x**i for i in range(fit_polynom_order + 1)]))


odr_model = odr.Model(func)
odr_data = odr.Data(time_x, time_y)
odr_obj = odr.ODR(odr_data, odr_model, beta0=np.zeros(fit_polynom_order + 1))
odr_res = odr_obj.run()
par_best = odr_res.beta

# Symbolical representation of the function (used by the sympy module):
function = np.dot(par_best, [x**i for i in range(len(par_best))])


""" ---------------------------------------------------------------------- """

"""
    Calculate the projections of the datapoints on the curve
"""


def r(x, y, Xi, Yi):
    """ Return the distance between (x, y) and (Xi, Yi) """
    return ((x-Xi)**2 + (y - Yi)**2)**0.5


l_function = sp.lambdify(x, function)

def point_on_curve(Xi, Yi):
    """
        Return the (x, y)-coordinates of the point on the curve closest
        to (Xi, Yi)
    """
    global l_function

    # Pythagoras's theorem gives us the distance:
    distance = sp.real_root((Xi - x)**2 + (Yi - function)**2, 2)
    # We're interested in the points where the derivative of the distance equals zero:
    derivative = sp.fraction(distance.diff().simplify())
    derivative_zero = sp.solveset(derivative[0], x)
    # Previous line returns all solutions in the complex plane while we only want the real answers:
    derivative_zero = np.array([n if sp.re(n) == n else None
                                for n in derivative_zero])
    derivative_zero = derivative_zero[derivative_zero != None]

    # Return the result closest to the function
    shortest_distance = np.Inf
    x_new = np.Inf
    for x_i in derivative_zero:
        if r(x_i, l_function(x_i), Xi, Yi) < shortest_distance:
            shortest_distance = r(x_i, l_function(x_i), Xi, Yi)
            x_new = x_i

    return (x_new, l_function(x_new))


# make an array with the projections of the datapoints on the curve
data_on_curve = np.array([point_on_curve(data_x[i], data_y[i])
                          for i in range(len(data_x))]).T


""" ---------------------------------------------------------------------- """

'''
    Write data to output file (lines beginning with a '#' are comments)
'''

output_file = open(output_filename, 'w')
output_file.write('# Data projected on timecurve \n')
output_file.write('# Timecurve fitted to the datapoints from the time series\n')
output_file.write('# Timecurve: {}\n'.format(function))                
output_file.write('# Gene \t\t\t\t PC1 \t\t PC2 \t x on curve\t\ty on curve \n')
for i in range(len(data_x)):
    output_file.write('{} \t\t {:10.2f} \t {:7.2f} \t {:7.2f} \t {:10.2f} \n'.format(gene[i].decode('UTF-8'), data_x[i], data_y[i], data_on_curve[0][i], data_on_curve[1][i]))
    
output_file.close()


""" ---------------------------------------------------------------------- """

''' Plot the results '''
if savefig:
    ax = plt.subplot(111)
    x_curve = np.linspace(min(time_x) - 2, max(time_x) + 2, 1000)
    ax.plot(x_curve, func(par_best, x_curve), label="Timecurve")
    ax.plot(data_x, data_y, '.', label="Original datapoints", linewidth=0.5, color='crimson')
    ax.plot(time_x, time_y, '.', label="Timeseries", color="navy")
    ax.plot(data_on_curve[0], data_on_curve[1], '.', label="Projection on curve", linewidth=0.5, color='green')
    plt.legend()
    plt.title("Projection on timecurve")
    ax.set_aspect('equal')
    plt.tight_layout()
    plt.savefig(inputfile_figure, format='pdf')
