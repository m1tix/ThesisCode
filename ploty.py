#!/usr/bin/python3
import csv
import sys
import json

import matplotlib.pyplot as plt
import numpy as np


def read_dat(dat_file):
    x_axes = []
    algebraic_method = []
    fincke_method = []
    with open(dat_file, "r") as f:
        normaldat = csv.reader(f)
        for row in normaldat:
            x_axes.append(int(row[0]))
            algebraic_method.append(float(row[1]))
            fincke_method.append(float(row[2]))
    return x_axes, algebraic_method, fincke_method


def plotly():
    x_axes, algebraic_method, fincke_method = read_dat("data/dat1")
    plt.figure()
    plt.scatter(x_axes,
                algebraic_method,
                color='red',
                s=15,
                label="Class group method")
    plt.scatter(x_axes,
                fincke_method,
                color='blue',
                s=15,
                label="Fincke-Pohst")
    plt.title(
        r"Time (s) to determine $|N_{\mathbb{Q}(\sqrt{506})/\mathbb{Q}}(x)| = m$ for all integral $x$."
    )
    plt.xlabel('m')
    plt.ylabel('time (s)')
    plt.legend(loc='upper left')
    plt.show()


def plotly2():
    x_axes, algebraic_method, fincke_method = read_dat("data/dat2")
    plt.figure()
    plt.plot(np.log10(x_axes),
             algebraic_method,
             color='red',
             label="Class group method")
    plt.plot(np.log10(x_axes),
             fincke_method,
             color='blue',
             label="Fincke-Pohst")
    plt.xlabel('n')
    plt.ylabel('time (s)')
    plt.title(
        r"Time (s) to calculate $|N_{\mathbb{Q}(\sqrt{506})/\mathbb{Q}}(x)| \leq \lfloor10^n\rfloor$ for all integral $x$."
    )
    plt.legend(loc='upper left')
    plt.show()

def plotly3():
    x_axes1, algebraic_method1, fincke_method1 = read_dat("data/dat1")
    with open("data/bak.json", "r") as f:
        raw_data = json.load(f)
    x_axes2 = []
    algebraic_method2 = []
    fincke_method2 = []
    for data_point in raw_data:
        x_axes2.append(data_point[0])
        algebraic_method2.append(np.mean(data_point[1]))
        fincke_method2.append(np.mean(data_point[2]))
    fig, (ax1, ax2) = plt.subplots(1, 2)
    fig.suptitle(
        r"Comparison of the two methods for $L=\mathbb{Q}(\sqrt{506})$",
        size=11, fontweight='bold')
    ax1.scatter(x_axes1,
                algebraic_method1,
                color='red',
                s=5,
                label="Class group method")
    ax1.scatter(x_axes1,
                fincke_method1,
                color='blue',
                s=5,
                label="Fincke-Pohst method")
    ax2.scatter(np.log10(x_axes2),
             algebraic_method2,
             color='red',
             s=5,
             label="Class group method")
    ax2.scatter(np.log10(x_axes2),
             fincke_method2,
             color='blue',
             s=5,
             label="Fincke-Pohst method")
    ax1.set_title(r"Time (s) to calculate $|N_{L/\mathbb{Q}}(x)| = m$", size=9)
    ax2.set_title(
        r"Time (s) to calculate $|N_{L/\mathbb{Q}}(x)| \leq \lfloor10^n\rfloor$",
        size=9)
    ax1.legend(loc='upper left', fontsize=6)
    ax2.legend(loc='upper left', fontsize=6)
    ax1.set(xlabel=r'$m$', ylabel=r'time (s)')
    ax2.set(xlabel=r'$n$', ylabel=r'time (s)')
    plt.tight_layout()
    plt.show()

def plotly4():
    with open("data/bak.json", "r") as f:
        raw_data = json.load(f)
    x_axes = []
    y_axes = []
    y_err = []
    for data_point in raw_data:
        x_axes.append(data_point[0])
        y_axes.append(np.mean(data_point[1]))
        y_err.append(np.std(data_point[1]))
    plt.errorbar(np.log10(x_axes), y_axes, yerr=y_err)
    plt.show()

if __name__ == "__main__":
    args = sys.argv
    if len(args) == 1:
        print("Specify what to plot by passing 1, 2 or 3")
        sys.exit(0)
    else:
        if args[1] == '1':
            plotly()
        elif args[1] == '2':
            plotly2()
        elif args[1] == '3':
            plotly3()
        else:
            print("Unreadable input")
            sys.exit(0)
