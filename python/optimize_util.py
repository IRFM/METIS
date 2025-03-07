#!/usr/bin/env python

import numpy as np
import pandas as pd
import seaborn
import os
import matplotlib as mpl
import matplotlib.pyplot as plt
from bayes_opt import BayesianOptimization
from bayes_opt.bayesian_optimization import acq_max, unique_rows
import datetime
import warnings
from collections import OrderedDict
from pymetis import matfunc2py, import_function


MAX_COST_VALUE=10


def xy_df_plot(df,
               xs,
               ys,
               ncols=4,
               xlim=(None, None),
               ylim=(None, None),
               **plotkwargs):

    import matplotlib.pyplot as plt
    import itertools

    if isinstance(xs, str):
        xs = [xs]
    if isinstance(ys, str):
        ys = [ys]
    n = len(xs) * len(ys)
    nr = n // ncols + int((n % ncols) > 0)
    fig, axs = plt.subplots(nr, min(ncols, n), **plotkwargs)
    if n == 1:
        axs = np.array([axs])
    for ((x, y), ax) in zip(itertools.product(xs, ys), axs.flatten()):
        plt.sca(ax)
        df.plot(kind='scatter', x=x, y=y, ax=ax)
        limits = list(ax.get_xlim())
        if xlim[0] is not None:
            limits[0] = xlim[0]
        if xlim[1] is not None:
            limits[1] = xlim[1]
        else:
            limits[1] = df[x].max()
        ax.set_xlim(limits)
        limits = list(ax.get_ylim())
        if ylim[0] is not None:
            limits[0] = ylim[0]
        if ylim[1] is not None:
            limits[1] = ylim[1]
        ax.set_ylim(limits)


def plot_0d(df,
            metisx,
            fields=(('cons', ('ip', 'pecrh', 'xece')),
                    ('zerod',
                     ('li', 'q0', 'q95', 'te0', 'pin', 'prad', 'plhthr'))),
            posts_all=None,
            posts_spec=None,
            metis_fast=False,
            it1=1,
            legend=False,
            titles=False,
            fname=None):

    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import itertools

    y_ref = None
    if posts_spec is not None:
        for label, post in posts_spec:
            if label == 'reference':
                y_ref = (('nbar', post.z0dinput.cons.nbar * 1e-19),
                         ('pnbi', (np.real(post.z0dinput.cons.pnbi) +
                                   np.imag(post.z0dinput.cons.pnbi)) * 1e-6),
                         ('picrh', np.abs(post.z0dinput.cons.picrh * 1e-6)), )
                post_ref = post

    for grp, flds in fields:

        nrows = len(flds) // 2 + len(flds) % 2
        if y_ref is None or grp != 'cons':
            fig, axs = plt.subplots(nrows=nrows,
                                    ncols=2,
                                    sharex=True,
                                    figsize=(15, 10))
        else:
            nrows += (len(flds) + 1) % 2
            ly = len(y_ref)
            fig = plt.figure(figsize=(15, 10))
            gs = gridspec.GridSpec(ly * nrows, 2)
            axs = []
            for i, ii in zip(range(len(flds) + 1), itertools.cycle((0, 1))):
                if i < len(flds):
                    axs.append(plt.subplot(gs[(i // 2) * ly:((i // 2) + 1) *
                                              ly, ii:ii + 1]))
                    axs[-1].tick_params('x', labelbottom=False)
                    axs[-1].locator_params('y', nbins=4)
                else:
                    # for iy in range(len(y_ref)):
                    #     axs.append(plt.subplot(gs[(i//2)*ly+iy:(i//2)*ly+iy+1, ii:ii+1]))
                    for iy, y in enumerate(y_ref):
                        ax = plt.subplot(gs[(i // 2) * ly + iy:(i // 2) * ly +
                                            iy + 1, ii:ii + 1])
                        ax.plot(post_ref.z0dinput.cons.temps,
                                y[1],
                                label=y[0],
                                c='k')
                        if legend:
                            ax.legend(loc='upper left')
                        ylim = ax.get_ylim()
                        ax.tick_params('x', labelbottom=False)
                        ax.locator_params('y', nbins=2)
                        if ylim[1] <= y[1].max() * 1.05:
                            ax.set_ylim(ylim[0], y[1].max() * 1.05)
                    ax.tick_params('x', labelbottom=True)
            axs = np.array(axs)

        for di, d in df.iterrows():
            if posts_all is None:
                post = metisx(return_post=True,
                              fast=metis_fast,
                              **{k: v
                                 for k, v in d.iteritems()
                                 if k.startswith('x')})
            else:
                post = posts_all[di]
            if grp == 'cons':
                data = post['z0dinput']['cons']
            else:
                data = post[grp]

            for ax, fld in zip(axs.reshape(-1), flds):
                ydata = data[fld][it1:]
                if fld == 'te0':
                    ydata = ydata * 1e-3
                if fld in ('ip', 'pecrh'):
                    ydata = ydata * 1e-6
                ax.plot(data['temps'][it1:], ydata, label=str(di))
                if titles:
                    ax.set_title(fld)
                # limits
                if fld == 'li':
                    ax.axhline(0.5, ls='--', c='k')
                    ax.axhline(1.2, ls='--', c='k')
                if fld == 'qmin':
                    ax.axhline(1, ls='--', c='k')
                ax.tick_params('x', bottom=False)
                ax.locator_params('y', nbins=4)
            ax.tick_params('x', labelbottom=True)

        if posts_spec is not None:
            for label, post in posts_spec:
                if grp == 'cons':
                    data = post['z0dinput']['cons']
                else:
                    data = post[grp]
                for ax, fld in zip(axs.reshape(-1), flds):
                    ydata = data[fld][it1:]
                    if fld == 'te0':
                        ydata = ydata * 1e-3
                    if fld in ('ip', 'pecrh'):
                        ydata = ydata * 1e-6
                    ax.plot(data['temps'][it1:], ydata, '--', label=label)
                    if fld == 'flux_tot':
                        ylim = ax.get_ylim()
                        ax.set_ylim(ydata.max() - 80, ydata.max() + 5)
                    if fld == 'li':
                        ax.set_ylim(0, 2)
                    if fld == 'qmin':
                        ax.set_ylim(0, 3)
                    if fld == 'q95':
                        ax.set_ylim(0, 7)

                    ax.tick_params('x', labelbottom=False)
                    # ax.locator_params('y', nbins=4)
                ax.tick_params('x', labelbottom=True)

        if grp != 'cons':
            axs.reshape(-1)[-1].tick_params('x', labelbottom=True)
            axs.reshape(-1)[-2].tick_params('x', labelbottom=True)

        if legend:
            ax.legend(loc='best')

        if fname is not None:
            fname2 = os.path.splitext(fname)[0] + grp + os.path.splitext(
                fname)[1]
            fig.savefig(fname2, dpi=200, bbox_inches='tight')


def plot_1d(df,
            metisx,
            fields=('qjli', 'jli', 'tep', 'tip'),
            tidx=-1,
            posts_all=None,
            posts_spec=None,
            metis_fast=False,
            ncols=2,
            fname=None,
            **plotkwargs):

    import matplotlib.pyplot as plt

    n = len(fields)
    nr = n // ncols + int((n % ncols) > 0)
    fig, axs = plt.subplots(nr, min(ncols, n), sharex=True, **plotkwargs)

    if posts_spec is None:
        posts = []
        lstyles = []
    else:
        posts = list(posts_spec)
        lstyles = ['--' for post in posts]

    for di, d in df.iterrows():
        if posts_all is None:
            post = metisx(return_post=True,
                          fast=metis_fast,
                          **{k: v
                             for k, v in d.iteritems() if k.startswith('x')})
        else:
            post = posts_all[di]
        posts.insert(0, (str(di), post))
        lstyles.insert(0, '-')

    for fld, ax in zip(fields, axs.reshape(-1)):
        for ((label, post), ls) in zip(posts, lstyles):
            ydata = post['profil0d'][fld][tidx]
            if fld in ('tip', 'tep'):
                ydata = ydata * 1e-3
            if fld in ('jli', ):
                ydata = ydata * 1e-6

            ax.plot(post['profil0d']['xli'][0], ydata, ls=ls, label=label)

        h, l = ax.get_legend_handles_labels()

        # ax.set_title(fld)

    figl = plt.figure(figsize=(3, 9))
    plt.figlegend(h, l, "upper right")

    if fname is not None:
        fig.savefig(fname, dpi=200, bbox_inches='tight')
        fname2 = os.path.splitext(fname)[0] + '_legend' + os.path.splitext(
            fname)[1]
        figl.savefig(fname2, dpi=200, bbox_inches='tight')

    # ax.legend(loc='best')


def bo_res_df(bo_res, posts_all, run=1):
    # construct DataFrame from posts_all
    df = pd.DataFrame(bo_res['all']['params'])

    df['J_bo'] = bo_res['all']['values']
    df['run'] = run

    # calculate all the objective funtions
    df['J_ss'] = [JD_ss(post, 0) for post in posts_all]
    df['J_sq'] = [JD_sq(post, 0) for post in posts_all]
    df['J_co'] = [JD_combined(post, 0) for post in posts_all]

    df['cons_D'] = [J_constraints_Dongen(post) for post in posts_all]
    df['cons_li'] = [J_constraints_Dongen(post,
                                          q0=False,
                                          ramp=False,
                                          li=True,
                                          hmode=False) for post in posts_all]
    df['cons_q0'] = [J_constraints_Dongen(post,
                                          q0=True,
                                          ramp=False,
                                          li=False,
                                          hmode=False) for post in posts_all]
    df['cons_ramp'] = [J_constraints_Dongen(post,
                                            q0=False,
                                            ramp=True,
                                            li=False,
                                            hmode=False) for post in posts_all]
    df['cons_hmode'] = [J_constraints_Dongen(post,
                                             q0=False,
                                             ramp=False,
                                             li=False,
                                             hmode=True) for post in posts_all]

    df['J_Dongen_co'] = [JD_combined(post) for post in posts_all]
    df['J_Dongen_ss'] = [JD_ss(post) for post in posts_all]
    df['J_Dongen_sq'] = [JD_sq(post) for post in posts_all]

    return df


def posterior2(gp, X, Y, x):
    ur = unique_rows(X)
    gp.fit(X[ur], Y[ur])
    # gp.fit(X, Y)
    try:
        mu, sigma2 = gp.predict(x.reshape(-1, X[0].size), eval_MSE=True)
        sigma = np.sqrt(sigma2)
        warnings.warn('detected obsolet GaussianProcess')
    except Exception:
        mu, sigma = gp.predict(x.reshape(-1, X[0].size), return_std=True)

    return mu, sigma


def plot_gp(bo, x, y=None, resfig=False, niter=None, ylim=(None, None),
            xlabel='x',
            ylabel1='Objective function',
            ylabel2='Utility'):

    X, Y = bo.X[:niter], bo.Y[:niter]

    fig = plt.figure(figsize=(16, 10))
    # fig.suptitle('Optimization iteration {}'.format(len(X)))
        # fontdict={'size': 50})

    gs = mpl.gridspec.GridSpec(2, 1, height_ratios=[3, 1])
    axis = plt.subplot(gs[0])
    axis.set_title('Optimization iteration {}'.format(len(X)))
    acq = plt.subplot(gs[1])

    mu, sigma = posterior2(bo.gp, X, Y, x)
    if y:
        axis.plot(x, y, linewidth=3, label='Target')
    axis.plot(X[:-1].flatten(),
              Y[:-1],
              'D',
              markersize=8,
              label=u'Observations',
              color='r')
    axis.plot(X[-1].flatten(),
              Y[-1],
              'x',
              markersize=15,
              markerfacecolor='k',
              markeredgecolor='k',
              markeredgewidth=1,
              label=u'Last samle')
    axis.axhline(Y.max(), ls=':', lw=2, color='b', label='Maximum')
    axis.plot(x, mu, '--', color='k', label='Prediction')

    axis.fill(
        np.concatenate([x, x[::-1]]),
        np.concatenate([mu - 1.9600 * sigma, (mu + 1.9600 * sigma)[::-1]]),
        alpha=.6,
        fc='c',
        ec='None',
        label='95% confidence interval')

    axis.set_xlim(x[0], x[-1])
    axis.set_ylim(ylim)
    axis.set_ylabel(ylabel1)  # , fontdict={'size': 30})
    axis.set_xlabel(xlabel)  # , fontdict={'size': 30})

    utility = bo.util.utility(x.reshape((-1, 1)), bo.gp, Y.max())

    acq.plot(x, utility, label='Utility Function', color='purple')

    # Update maximum value to search for next probe point.
    y_max = Y.max()

    # Maximize acquisition function to find next probing point
    if niter < bo.X.shape[0]:
        # take x_max from the next iteration
        x_max = bo.X[niter]
    else:
        x_max = acq_max(ac=bo.util.utility,
                        gp=bo.gp,
                        y_max=y_max,
                        bounds=bo.bounds)
    # acq.plot(x[np.argmax(utility)],
    #          np.max(utility),
    axis.plot(x_max,
              y_max,
              '*',
              markersize=15,
              label=u'Next Best Guess',
              markerfacecolor='gold',
              markeredgecolor='k',
              markeredgewidth=1)
    acq.plot(x_max,
             bo.util.utility(x_max.reshape(-1, 1), bo.gp, Y.max()),
             '*',
             markersize=15,
             label=u'Next Best Guess',
             markerfacecolor='gold',
             markeredgecolor='k',
             markeredgewidth=1)
    axis.set_xlim(x[0], x[-1])
    #     acq.set_ylim((0, np.max(utility) * 1.1))
    acq.set_ylabel(ylabel2)  # , fontdict={'size': 30})
    acq.set_xlabel(xlabel)  # , fontdict={'size': 30})

    axis.legend(loc=2, bbox_to_anchor=(1.01, 1), borderaxespad=0.)
    acq.legend(loc=2, bbox_to_anchor=(1.01, 1), borderaxespad=0.)

    if resfig:
        return fig


def plot_gp_evolution(bo,
                      figdir,
                      nmin=None,
                      nmax=None,
                      x_plot=None,
                      ext="png",
                      timestamp=True,
                      ylim=(None, None),
                      xlabel='x',
                      ylabel1='Objective function',
                      ylabel2='Utility'):
    if timestamp is True:
        timestamp = '{:%Y-%m-%d-%H-%M}'.format(datetime.datetime.now())
    if not timestamp:
        timestamp = ''

    if x_plot is None:
        x_plot = np.linspace(bo.X.min(), bo.X.max(), 201)

    if nmin is None:
        nmin = len(bo.init_points)

    if nmax is None:
        nmax = bo.X.shape[0]
    nmax = min(nmax, bo.X.shape[0])

    for niter in range(nmin, nmax):
        name = 'bo_evol_{timestamp}_{niter:04d}.{ext}'.format(
            timestamp=timestamp,
            niter=niter,
            ext=ext)
        filename = os.path.join(figdir, name)
        fig = plot_gp(bo, x_plot, resfig=True, niter=niter, ylim=ylim,
                      xlabel=xlabel, ylabel1=ylabel1, ylabel2=ylabel2)
        fig.savefig(filename, bbox_inches='tight', pad_inches=0.1)


def cartesian(arrays, out=None):
    """
    Generate a cartesian product of input arrays.

    Parameters
    ----------
    arrays : list of array-like
        1-D arrays to form the cartesian product of.
    out : ndarray
        Array to place the cartesian product in.

    Returns
    -------
    out : ndarray
        2-D array of shape (M, len(arrays)) containing cartesian products
        formed of input arrays.

    Examples
    --------
    >>> cartesian(([1, 2, 3], [4, 5], [6, 7]))
    array([[1, 4, 6],
           [1, 4, 7],
           [1, 5, 6],
           [1, 5, 7],
           [2, 4, 6],
           [2, 4, 7],
           [2, 5, 6],
           [2, 5, 7],
           [3, 4, 6],
           [3, 4, 7],
           [3, 5, 6],
           [3, 5, 7]])

    """

    arrays = [np.asarray(x) for x in arrays]
    dtype = arrays[0].dtype

    n = np.prod([x.size for x in arrays])
    if out is None:
        out = np.zeros([n, len(arrays)], dtype=dtype)

    m = n // arrays[0].size
    out[:, 0] = np.repeat(arrays[0], m)
    if arrays[1:]:
        cartesian(arrays[1:], out=out[0:m, 1:])
        for j in range(1, arrays[0].size):
            out[j*m:(j+1)*m, 1:] = out[0:m, 1:]
    return out


def plot_gp_2D(bo, nx=51, ny=41, resfig=False, niter=None, ylim=(None, None),
               xlabel='', ylabel='', nc=5, cfmt='%.2f', legend_loc='upper right', X=None, Y=None):

    if X is None and Y is None:
        X_in = bo.X
        X, Y = bo.X[:niter], bo.Y[:niter]

        x1 = np.linspace(bo.X[:, 0].min(), bo.X[:, 0].max(), nx)
        x2 = np.linspace(bo.X[:, 1].min(), bo.X[:, 1].max(), ny)
    else:
        X_in = X
        # take limits from all iterations
        x1 = np.linspace(X[:, 0].min(), X[:, 0].max(), nx)
        x2 = np.linspace(X[:, 1].min(), X[:, 1].max(), ny)
        # select the iteration to plot
        X, Y = X[:niter], Y[:niter]

    xx = cartesian((x1.squeeze(), x2.squeeze()))
    X1, X2 = np.meshgrid(x1, x2)

    ur = unique_rows(X)
    bo.gp.fit(X[ur], Y[ur])
    # gp.fit(X, Y)
    mu, sigma = posterior2(bo.gp, X, Y, xx)
    utility = bo.util.utility(xx, bo.gp, Y.max())

    # Maximize acquisition function to find next probing point
    if niter is not None and niter < X_in.shape[0]:
        # take x_max from the next iteration
        x_max = X_in[niter]
    else:
        x_max = acq_max(ac=bo.util.utility,
                        gp=bo.gp,
                        y_max=Y.max(),
                        bounds=bo.bounds)

    # cmapf = plt.cm.RdPu_r
    cmapf = seaborn.light_palette("blue", n_colors=6, reverse=False, as_cmap=True)
    plot_2d(X1, X2, F1=mu, F2=sigma, X=X, Y=Y,
            x_max=x_max,
            title='Objective function prediction (lines) and standard deviation (color)',
            xlabel=xlabel,
            ylabel=ylabel,
            nc=nc,
            cmapf=cmapf,
            cfmt=cfmt,
            legend_loc=legend_loc)

    plot_2d(X1, X2, F1=mu, F2=utility, X=X, Y=Y,
            x_max=x_max,
            title='Objective function prediction (lines) and utility function (color)',
            xlabel=xlabel,
            ylabel=ylabel,
            nc=nc,
            cmapf=cmapf,
            cfmt=cfmt,
            legend_loc=legend_loc)

    print("x_max = {}".format(x_max))

    if resfig:
        return plt.gcf()


def plot_2d(X1, X2, F1, F2, X, Y,
            x_max=None,
            title='',
            xlabel='',
            ylabel='',
            nc=5,
            cmapf=plt.cm.RdPu_r,
            cfmt='%.2f',
            legend_loc='upper right'):
        plt.figure()
        # cmap = plt.cm.spring_r
        cmap = seaborn.dark_palette("red", n_colors=nc, reverse=False, as_cmap=True)
        c = plt.contour(X1, X2, F1.reshape(X1.shape[::-1]).T, nc,
                        cmap=cmap, label='Target')
        # plt.clabel(c, inline=1, fontsize=10)

        plt.clabel(c, fmt=cfmt)
    #     plt.colorbar(c)
        plt.scatter(X[:, 0], X[:, 1], c=Y, cmap=cmap)

    #     plt.figure()
        c = plt.contourf(X1, X2, F2.reshape(X1.shape[::-1]).T, nc,
                         cmap=cmapf,
                         alpha=0.3, label='Standard deviation')
        plt.colorbar(c)
        imax = Y.argmax()
        plt.plot(X[-1, 0], X[-1, 1],
                 ls='',
                 marker='x',
                 markersize=10,
                 label=u'Last Sample',
                 markerfacecolor='k',
                 markeredgecolor='k',
                 markeredgewidth=2)
        plt.plot(X[imax, 0], X[imax, 1],
                 ls='',
                 marker='+',
                 markersize=10,
                 label=u'Maximum',
                 markerfacecolor='blue',
                 markeredgecolor='b',
                 markeredgewidth=2)
        if x_max is not None:
            plt.plot(x_max[0],
                     x_max[1],
                     ls='',
                     marker='*',
                     markersize=10,
                     label=u'Next Best Guess',
                     markerfacecolor='g',
                     markeredgecolor='g',
                     markeredgewidth=2)
        if legend_loc:
            plt.legend(loc=legend_loc, shadow=True, fancybox=True, frameon=True)
        plt.title(title)
        plt.xlabel(xlabel)
        plt.ylabel(ylabel)


def initialize_from_res(res, metisx, normalize=False):
    bo = BayesianOptimization(metisx, metisx.get_x_limits())
    # transform the res dict to the required form
    ymin = min(res['all']['values'])
    ymax = max(res['all']['values'])
    print('ymin = {}, ymax = {}'.format(ymin, ymax))
    if normalize:
        metisx.normalize(ymin, ymax)

        def fnorm(cost):
            return norm_cost(cost, ymin, ymax)
    else:

        def fnorm(cost):
            return cost

    points_dict = {fnorm(y): x
                   for y, x in zip(res['all']['values'], res['all']['params'])}

    bo.initialize(points_dict)

    return bo


def new_objective(self, obj_func):
    """Use a new objective function, recalculate results
    """
    self.f = obj_func
    Y = [self.f(**dict(zip(self.keys, x))) for x in self.X]
    self.Y = np.array(Y)
    self.res['all']['params'] = [dict(zip(self.keys, x))
                                 for x in self.X[len(self.init_points):]]
    self.res['all']['values'] = [y for y in self.Y[len(self.init_points):]]
    imax = self.Y.argmax()
    self.res['max']['max_params'] = dict(zip(self.keys, self.X[imax]))
    self.res['max']['max_val'] = self.Y[imax]


def explore_new(self, X):
        # X is a numpy array
    #     i0 = self.X.shape[0]
    #     self.X = np.resize(self.X, (i0 + len(points_dict), self.X.shape[1]))
    #     self.Y = np.resize(self.Y, (self.X.shape[0], ))
    #     for p, Xi, Yi in zip(points_dict, self.X[i0: ]):
    #         Xi[()] = np.array([p[k] for k in self.keys])
    #         Yi[()] = self.f(**points_dict)

    #         self.res['all']['params'].append(
    #     self.res['all']['values'] = [y for y in self.Y[len(self.init_points): ]]
    #     imax = self.Y.argmax()
    #     self.res['max']['max_params'] = dict(zip(self.keys, self.X[imax]))
    #     self.res['max']['max_val'] = self.Y[imax]

        self.X = np.vstack((self.X, X))
        Y_new = [self.f(**dict(zip(self.keys, x))) for x in X]
        self.Y = np.hstack((self.Y, Y_new))
        for x, y in zip(X, Y_new):
            self.res['all']['params'].append(dict(zip(self.keys, x)))
            self.res['all']['values'].append(y)
        imax = self.Y.argmax()
        self.res['max']['max_params'] = dict(zip(self.keys, self.X[imax]))
        self.res['max']['max_val'] = self.Y[imax]


# patch the original class
BayesianOptimization.new_objective = new_objective
BayesianOptimization.explore_new = explore_new


class NormFunc(object):
    """Automatically normalized function
    """

    def __init__(self, func, fmax=1, fmin=0):
        self.fmax = fmax
        self.fmin = fmin
        self.func = func

    def __call__(self, *args, normalize=True, **kwargs):
        res = self.func(*args, **kwargs)

        if normalize:
            res = (res - self.fmin) / (self.fmax - self.fmin)

        return res

    def revert(self, value):
        """Convert normalized value to un-normalized
        """
        res = value * (self.fmax - self.fmin) + self.fmin
        return res


class Constraint(object):

    def __init__(self, func, fmin=None, fmax=None, power=1, mult=1, ulim=None):
        assert fmin is not None or fmax is not None
        self.func = func
        self.fmin = fmin
        self.fmax = fmax
        self.power = power
        self.mult = mult
        self.ulim = ulim

    def __call__(self, *args, return_value=False, **kwargs):
        value = self.func(*args, **kwargs)
        if self.fmin is not None:
            overmin = np.maximum(self.fmin - value, 0)
        else:
            overmin = 0
        if self.fmax is not None:
            overmax = np.maximum(value - self.fmax, 0)
        else:
            overmax = 0

        penalty = (self.mult * (overmin + overmax))**self.power

        if self.ulim is not None:
            penalty = np.minimum(penalty, self.ulim)

        if return_value:
            return penalty, value
        else:
            return penalty


class MultiObjective(object):

    def __init__(self, objectives=(), constraints=()):
        self.obj_funcs = OrderedDict(objectives)
        self.cons_funcs = OrderedDict(constraints)

    def __call__(self, *args, eval_obj=True, eval_cons=True, **kwargs):
        assert eval_obj or eval_cons
        res = 0
        if eval_obj:
            for f, w in self.obj_funcs.items():
                res += w * f(*args, **kwargs)
        if eval_cons:
            for f, w in self.cons_funcs.items():
                res -= w * f(*args, **kwargs)

        return res


def norm_cost(cost, ymin=0, ymax=1):
    return (cost - ymin) / (ymax - ymin)


def J_Dongen(post, nu_ss=0, factor=1):
    '''Hybrid mode cost function

    J. van Dongen, F. Felici, G. M. D. Hogeweij, P. Geelen, and E. Maljaars,
    Plasma Phys. Control. Fusion 56, 125008 (2014). formula (4)
    '''
    from scipy.interpolate import UnivariateSpline
    from scipy.integrate import cumtrapz, trapz

    # s/q factor for hybrid modes [J. Citrin et al., Plasma Phys. Control. Fusion 54 (2012).]
    # (2) J_sq = -<s/q>
    rho = np.sqrt(post.profil0d.phi[-1, :])
    rho /= rho[-1]  # normalize
    # q_rho = interp1d(rho, post.profil0d.qjli[-1,:], kind=3)
    q_rho = UnivariateSpline(rho, post.profil0d.qjli[-1, :])
    q = q_rho(rho)
    dqdrho = q_rho.derivative()(rho)
    s = rho / q * dqdrho
    V = cumtrapz(post.profil0d.vpr[-1, :], post.profil0d.xli[-1, :], initial=0)
    V_rho = UnivariateSpline(rho, V)
    dVdrho = V_rho.derivative()(rho)

    Jsq = -trapz(dVdrho * s / q, rho) / V[-1]

    # steady-state
    # (3) J_ss
    Wss = (1 - rho)
    Upar = post.profil0d.epar[-1, :] * (2 * np.pi * post.z0dinput.geo.R[-1])
    Jss = trapz((Wss * (Upar - Upar[-1]))**2, rho)

    Jval = (1 - nu_ss) * Jsq + nu_ss * Jss
    Jval *= factor

    return Jval


def limit_cost_func(x,
                    limit,
                    power=10,
                    fact=1,
                    maxval=None,
                    absval=False,
                    type='semihard'):
    if absval:
        x = np.abs(x)
    if type == 'soft':
        res = fact * (np.exp((x - limit) * power))
    elif type == 'semihard':
        xa = np.asanyarray(x)
        res = np.zeros_like(xa)
        res[xa > limit] = fact * (np.exp((xa[xa > limit] - limit) * power) - 1)
        if res.size == 1 and not isinstance(x, np.ndarray):
            res = res.squeeze()[()]
    elif type == 'hard':
        res = np.zeros_like(x)
        res[x > limit] = np.infty
        if res.size == 1 and not isinstance(x, np.ndarray):
            res = res[()]
    elif type == 'sigmoid':
        res = maxval * 1 / (1 + np.exp(-(x - limit) * power))
    if maxval is not None:
        res = np.minimum(res, maxval)
    return res


def llimit_cost_func(x,
                     limit,
                     power=10,
                     fact=1,
                     maxval=None,
                     absval=False,
                     type='semihard'):
    return limit_cost_func(-x,
                           -limit,
                           power=power,
                           fact=fact,
                           maxval=maxval,
                           absval=absval,
                           type=type)


def fld_limit_cost_func(struct,
                        fld,
                        limit,
                        power=10,
                        fact=1,
                        maxval=None,
                        absval=False):
    return limit_cost_func(
        np.squeeze(struct[fld]),
        limit,
        power=power,
        fact=fact,
        maxval=maxval,
        absval=absval)


def ramp_up_limit_cost(z0dinput,
                       limit=0.2,
                       factor=0.5,
                       power=100,
                       maxval=100,
                       type='semihard'):
    # dIp / dt [MA/s]
    dIdt = np.diff(z0dinput.cons.ip.squeeze()) / np.diff(
        z0dinput.cons.temps.squeeze()) * 1e-6
    return factor * max(limit_cost_func(dIdt,
                                        limit,
                                        power=power,
                                        maxval=maxval,
                                        type=type))


def detect_error(func, retval=MAX_COST_VALUE):
    import functools

    @functools.wraps(func)
    def wrapper(post, *args, **kwargs):
        if 'error' in post:
            warnings.warn('error in post -> return nan')
            return retval
        else:
            return func(post, *args, **kwargs)

    return wrapper


@detect_error
def J_constraints_Dongen(post,
                         q0=True,
                         li=True,
                         ramp=False,
                         hmode=True,
                         fact=0.1):
    res = []

    # ramp-up limit
    if ramp:
        res.append(ramp_up_limit_cost(post.z0dinput, limit=0.2, power=30))

    if q0:
        # qmin > 1
        res.append(llimit_cost_func(post.zerod.qmin[3:].min(),
                                    1.0,
                                    fact=1,
                                    power=20))
        # q0 < 2
        res.append(0.5 * limit_cost_func(post.zerod.q0[-3:].mean(),
                                         2.0,
                                         power=20))

    if li:
        # li > 0.5
        # res.append(-llimit_cost_func(post.zerod.li[-1], 0.5, fact=2, power=20))
        res.append(2 * llimit_cost_func(post.zerod.li[3:].min(),
                                        0.5,
                                        power=20))
        # li < 1.2
        res.append(limit_cost_func(post.zerod.li[-1], 1.2, power=20))

    # H-mode enforced
    if hmode:
        pl2h = post.zerod.plhthr / post.zerod.plossl2h
        # last 3 time slices are hmodes
        res.append(llimit_cost_func(pl2h[-3:].min(), 1.2, power=10))

    return fact * sum(res).squeeze()[()]


@detect_error
def JD_combined(post, constraints=1, maxval=MAX_COST_VALUE, nanmax=False):
    res = []

    cons_fact = float(constraints)

    # s/q
    res.append(J_Dongen(post, 0))
    # steady-state
    res.append(10 * J_Dongen(post, 1))

    # constraints
    #     if constraints:
    res.append(cons_fact * J_constraints_Dongen(post))

    # sum up the total
    tot = sum(res)

    # filter out unexpectedly large values
    if tot > maxval:
        warnings.warn('maxval reached: %g > %g' % (tot, maxval))
        if nanmax:
            tot = np.nan
        else:
            tot = maxval

    return tot


@detect_error
def JD_sq(post, constraints=1, maxval=MAX_COST_VALUE, nanmax=False):
    res = []
    cons_fact = float(constraints)

    # s/q
    res.append(J_Dongen(post, 0))

    # constraints
    #     if constraints:
    res.append(cons_fact * J_constraints_Dongen(post))

    # sum up the total
    tot = sum(res)

    # filter out unexpectedly large values
    if tot > maxval:
        warnings.warn('maxval reached: %g > %g' % (tot, maxval))
        if nanmax:
            tot = np.nan
        else:
            tot = maxval

    return tot


@detect_error
def JD_ss(post, constraints=1, maxval=MAX_COST_VALUE, nanmax=False):
    res = []
    cons_fact = float(constraints)

    # steady-state
    res.append(10 * J_Dongen(post, 1))

    # constraints
    #     if constraints:
    res.append(cons_fact * J_constraints_Dongen(post))

    # sum up the total
    tot = sum(res)

    # filter out unexpectedly large values
    if tot > maxval:
        warnings.warn('maxval reached: %g > %g' % (tot, maxval))
        if nanmax:
            tot = np.nan
        else:
            tot = maxval

    return tot


z0fluxconsum = matfunc2py('z0fluxconsum', bridge=None, nargout=4)


@detect_error
def J_flux(post, constraints=1, maxval=None, nanmax=False):
    res = []

    cons_fact = float(constraints)

    flux_tot, flux_cs, flux_plasma_tot, flux_res = z0fluxconsum(
        post.z0dinput, post.zerod, post.profil0d)
    tot = flux_tot.squeeze()[-1]

    # filter out unexpectedly large values
    if maxval is not None and tot > maxval:
        warnings.warn('maxval reached: %g > %g' % (tot, maxval))
        if nanmax:
            tot = np.nan
        else:
            tot = maxval

    return tot


def penalty(func, g_min, g_max, alpha=1, beta=10, gamma=1):

    def output(*x):
        val = func(*x)
        if val < g_min:
            # res = min(alpha * (g_min - val) ** gamma, beta)
            res = alpha * (np.exp((g_min - val) * gamma) - 1)
        elif val > g_max:
            # res = beta / (1 + np.exp(-(val - g_max) * gamma))
            # res = min(alpha * (val - g_max) ** gamma, beta)
            res = alpha * (np.exp((val - g_max) * gamma) - 1)
        else:
            res = 0
            # res = beta / (1 + np.exp(-(val - g_max) * gamma))

        if beta is not None:
            res = np.minimum(res, beta)
        return res

    return output
