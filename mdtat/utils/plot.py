#!/usr/bin/env python
# AUTHOR:       Shane Gordon
# CREATED:      2015-06-16 21:46:32

import matplotlib.pyplot as plt
import numpy as np


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
    log = 'Attempting to plot to {}'.format(output)
    f, ax = plt.subplots(2)
    count = np.arange(0, len(data))
    for row, i in zip(data, count):
        ax[0].set_xlabel(xlabel)
        ax[0].set_ylabel(ylabel)
        ax[0].plot(row, label=i)
        # ax[0].set_ylim(0,)
        ax[0].legend()
        ax[1].set_xlabel(xlabel)
        ax[1].set_ylabel('ACF')
        ax[1].semilogx(autocorr(row), label=i)
        ax[1].legend()
    plt.savefig(output, format='pdf')
    plt.close()
    return log
