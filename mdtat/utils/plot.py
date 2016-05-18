#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit


def init_plotting():
    plt.rcParams['figure.figsize'] = (8, 3)
    plt.rcParams['figure.autolayout'] = True


def autocorr(x):
    "Compute an autocorrelation with numpy"
    x = x - np.mean(x)
    result = np.correlate(x, x, mode='full')
    result = result[result.size//2:]
    return result / result[0]


def basic_plot(data, xlabel, ylabel, output):

    def func(x, a, tau, offset):
        return a * np.exp(-x/tau) + offset

    def fit_autocorr(x, a=1, tau=1):
        y = func(x, a, tau)
        popt, pcov = curve_fit(func, x, y)
        return popt, pcov

    log = 'Attempting to plot to {}'.format(output)
    f, ax = plt.subplots(2, sharex=True)
    count = np.arange(0, len(data))
    for row, i in zip(data, count):
        acf = autocorr(row)
        x = np.arange(len(acf))
        popt, pcov = curve_fit(func, x, acf, p0=[1, 10, 0])
        ax[0].set_ylabel(ylabel)
        ax[0].plot(row, label=i)
        ax[0].legend()
        ax[1].set_xlabel(xlabel)
        ax[1].set_ylabel('ACF')
        ax[1].plot(acf, linestyle=':', label='Raw: {i}'.format(i=i))
        ax[1].plot(func(np.arange(len(acf)), popt[0], popt[1], popt[2]),
                   label='Fit: {rep}, tau={tau}, a={a}'.
                   format(rep=i, tau=popt[1], a=popt[0]))
        ax[1].legend(framealpha=0.75, loc=0, fontsize=6)
    plt.savefig(output, format='pdf')
    plt.close()
    return log
