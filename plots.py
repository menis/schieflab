#!/usr/bin/env python
# filename: plots.py


#
# Copyright (c) 2016 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#


import numpy as np
np.set_printoptions(suppress=True)
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


def shared_mutation_2dhist(x, y, cmap, ax, figfile=None, figsize=None, figsize_denom=4, n_levels=10, alpha=1.0, x_lim=(0, 14), y_lim=(0, 14),
                           y_label='VRC01-class mutations', x_label='Total amino acid mutations', labelsize=14,
                           tick_labelsize=12, pad=4, show_values=True):
    # adjust the inputs to make them iterable
    if type(x[0]) in [int, float, np.int64, np.float64]:
        x = [x]
        y = [y]
        cmap = [cmap]
    sns.set_style('whitegrid')
    # make the plots
    for _x, _y, _cmap in zip(x, y, cmap):
        bin_x = max(_x) + 2 if x_lim is None else x_lim[1] + 2
        bin_y = max(_y) + 2 if y_lim is None else y_lim[1] + 2
        bins = [[i - 0.5 for i in range(bin_y)], [i - 0.5 for i in range(bin_x)]]
        data, x_edges, y_edges = np.histogram2d(_y, _x, bins=bins)
        data = data[::1]
        if figsize is None:
            figsize = (float(bin_x) / figsize_denom, float(bin_y) / figsize_denom)
        mask = np.array([[val == 0 for val in subl] for subl in data])
        ax = sns.heatmap(data, cmap=_cmap, square=True, cbar=False, mask=mask,
                         linewidths=1, linecolor='w', alpha=alpha, annot=True, fmt=".0f")
    if x_lim is None:
        x_lim = [0, len(data[0])]
    else:
        x_lim = [x_lim[0], x_lim[1]]
    if y_lim is None:
        y_lim = [0, len(data)]
    else:
        y_lim = [y_lim[0], y_lim[1]]
    ax.set_xlim(x_lim)
    ax.set_ylim(y_lim)
    # format ticks and spines
    x_ticks = [t for t in range(x_lim[0], x_lim[1] + 1)]
    ax.set_xticks([t + 0.5 for t in x_ticks])
    x_ticklabels = [str(l) if l % 2 == 0 else '' for l in x_ticks]
    ax.set_xticklabels(x_ticklabels, rotation=0)
    y_ticks = [t for t in range(y_lim[0], y_lim[1] + 1)]
    ax.set_yticks([t + 0.5 for t in y_ticks])
    y_ticklabels = [str(l) if l % 2 == 0 else '' for l in y_ticks]
    ax.set_yticklabels(y_ticklabels, rotation=0)
    for position in ['right', 'left', 'top', 'bottom']:
        ax.spines[position].set_visible(True)
        ax.spines[position].set_color('k')
    if show_values:
        _data = data[::-1]
        max_val = float(np.amax(_data))
        for y in range(_data.shape[0]):
            for x in range(_data.shape[1]):
                if _data[y, x] == 0:
                    continue
                color = 'w' if _data[y, x] / max_val >= 0.5 else 'k'
                plt.text(x + 0.5, y + 0.5, str(int(_data[y, x])),
                         color=color,
                         horizontalalignment='center',
                         verticalalignment='center',
                        fontsize=12)
    ax.tick_params(axis='x', labelsize=tick_labelsize)
    ax.tick_params(axis='y', labelsize=tick_labelsize)
    plt.ylabel(y_label, size=labelsize)
    plt.xlabel(x_label, size=labelsize)
    ax.xaxis.set_tick_params(which='major', direction='out', length=3,
                             color='k', width=1, bottom=True, pad=pad)

    ax.yaxis.set_tick_params(which='major', direction='out', length=3,
                             color='k', width=1, left=True, pad=pad)


def pixel_plot(data, cmap, figfile=None, pad=2, labelsize=14, tight_layout=True,
               maxy_denom=30, maxx_denom=10):
    '''
    ::data:: is the output from vrc01_class_mutation_positions()
    '''
    max_y, max_x = data.shape
    f, ax = plt.subplots(figsize=(max_x / float(maxx_denom), max_y / float(maxy_denom)))
    plt.pcolor(data, cmap=cmap)
    ax.set_xlim([0, max_x])
    ax.set_ylim([0, max_y])
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_color('k')
    minorLocator = AutoMinorLocator(1)
    ax.xaxis.set_minor_locator(minorLocator)
    ticks = [26, 34, 50, 58, 99, 114]
    minor_ticks = [13, 30, 42, 54, 78.5, 106.5, 118]
    minor_labels = ['FR1', 'CDR1', 'FR2', 'CDR2', 'FR3', 'CDR3', 'FR4']
    ax.set_xticks(ticks)
    ax.xaxis.set_tick_params(which='major', direction='out',
                             length=6, color='k', width=1, top='off', labelsize=0)
    ax.set_xticks(minor_ticks, minor=True)
    ax.set_xticklabels(minor_labels, minor=True, y=-0.02)
    ax.xaxis.set_tick_params(which='minor', direction='out',
                             length=12, color='white', width=1, top='off',
                             labelsize=labelsize, pad=pad)
    ticks = [26, 34, 50, 58, 99, 114]
    plt.xticks(ticks, [' '] * len(ticks))
    ax.xaxis.set_tick_params(which='major', direction='out',
                             length=6, color='k', width=1.5, top='off')
    ax.tick_params(axis='y', labelsize=0)
    if tight_layout:
        plt.tight_layout()
    if figfile is None:
        plt.show()
    else:
        plt.savefig(figfile)


def fill_between_steps(ax, x, y1, y2, facecolor='grey', edgecolor='grey', alpha=0.5, step_where='mid', lw=1):
    ''' fill between a step plot.

    Parameters
    ----------
    ax : Axes
       The axes to draw to

    x : array-like
        Array/vector of index values.

    y1 : array-like or float
        Array/vector of values to be filled under.
    y2 : array-Like or float, optional
        Array/vector or bottom values for filled area. Default is 0.

    step_where : {'pre', 'post', 'mid'}
        where the step happens, same meanings as for `step`

    **kwargs will be passed to the matplotlib fill_between() function.

    Returns
    -------
    ret : PolyCollection
       The added artist

    '''
    if step_where not in {'pre', 'post', 'mid'}:
        raise ValueError("where must be one of {{'pre', 'post', 'mid'}} "
                         "You passed in {wh}".format(wh=step_where))
    # make sure y values are up-converted to arrays
    if np.isscalar(y1):
        y1 = np.ones_like(x) * y1
    if np.isscalar(y2):
        y2 = np.ones_like(x) * y2
    # temporary array for up-converting the values to step corners
    # 3 x 2N - 1 array
    vertices = np.vstack((x, y1, y2))
    if step_where == 'pre':
        steps = np.zeros((3, 2 * len(x) - 1), np.float)
        steps[0, 0::2], steps[0, 1::2] = vertices[0, :], vertices[0, :-1]
        steps[1:, 0::2], steps[1:, 1:-1:2] = vertices[1:, :], vertices[1:, 1:]
    elif step_where == 'post':
        steps = np.zeros((3, 2 * len(x) - 1), np.float)
        steps[0, ::2], steps[0, 1:-1:2] = vertices[0, :], vertices[0, 1:]
        steps[1:, 0::2], steps[1:, 1::2] = vertices[1:, :], vertices[1:, :-1]
    elif step_where == 'mid':
        steps = np.zeros((3, 2 * len(x)), np.float)
        steps[0, 1:-1:2] = 0.5 * (vertices[0, :-1] + vertices[0, 1:])
        steps[0, 2::2] = 0.5 * (vertices[0, :-1] + vertices[0, 1:])
        steps[0, 0] = vertices[0, 0]
        steps[0, -1] = vertices[0, -1]
        steps[1:, 0::2], steps[1:, 1::2] = vertices[1:, :], vertices[1:, :]
    else:
        raise RuntimeError("should never hit end of if-elif block for validated input")
    # un-pack
    xx, yy1, yy2 = steps
    # now to the plotting part:
    return ax.fill_between(xx, yy1, y2=yy2, facecolor='grey', edgecolor='grey', alpha=0.5, lw=1)




