""" Custom graphics lib for pop gen or genomics

FOREST Thomas (thomas.forest@college-de-france.fr)



"""

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import gc
import time
import datetime
import pandas as pd



def heatmap(data, row_labels=None, col_labels=None, ax=None,
            cbar_kw={}, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.
    (from the matplotlib doc)

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if not ax:
        ax = plt.gca()

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    if col_labels:
        ax.set_xticks(col_labels)
    if row_labels:
        ax.set_yticks(row_labels)
    
    # Let the horizontal axes labeling appear on top.
    ax.tick_params(top=True, bottom=False,
                   labeltop=True, labelbottom=False)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
             rotation_mode="anchor")

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1]+1)-.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0]+1)-.5, minor=True)
    ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar

def annotate_heatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
    """
    A function to annotate a heatmap.
     (from the matplotlib doc)
    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max())/2.

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center",
              verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
            text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
            texts.append(text)

    return texts

def plot_matrix(mat, legend=None, color_scale_type="YlGn", cbarlabel = "qt", title=None):
         
    fig, ax = plt.subplots(figsize=(10,8))
    if legend:
        row_labels = [k for k in range(len(mat))]
        col_labels = [k for k in range(len(mat[0]))]
        im, cbar = heatmap(mat, row_labels, col_labels, ax=ax,
                           cmap=color_scale_type, cbarlabel=cbarlabel)
    else:
        im, cbar = heatmap(mat, ax=ax,
                           cmap=color_scale_type, cbarlabel=cbarlabel)
    #texts = annotate_heatmap(im, valfmt="{x:.5f}")
    if title:
        ax.set_title(title)
    fig.tight_layout()
    plt.show()

def plot(x, y, outfile = None, outfolder = None, ylab=None, xlab=None,
         title=None, label = None, show=True, nb_subplots = None, subplot_init = False,
         subplot_id = None, output = None, dpi = 300, width = 15, height = 15, plot_init = True):

    # before fig is generated, set its dimensions
    if plot_init:
        plt.figure(figsize=(width, height))
    if subplot_init:
        # define a certain amount of subplots
        fig, axs = plt.subplots(nb_subplots)
    
    if x:
        if nb_subplots:
            axs[subplot_id].plot(x, y)
        else:
            fig, = plt.plot(x, y)
    else:
        # x is optional
        if nb_subplots:
            # define a certain amount of subplots
            axs[subplot_id].plot(y)
        else:
            fig, = plt.plot(y)
    if label:
        # if legend
        fig.set_label(label)
        plt.legend()
    if ylab:
        plt.ylabel(ylab)
    if xlab:
        plt.xlabel(xlab)
    if title:
        plt.title(title)
    if outfile:
        plt.savefig(outfile, dpi = dpi)
    if show == True:
        plt.show()


def scatter(x, y, ylab=None, xlab=None, title=None):
    plt.scatter(x, y)
    if ylab:
        plt.ylabel(ylab)
    if xlab:
        plt.xlabel(xlab)
    if title:
        plt.title(title)
    plt.show()

def barplot(x=None, y=None, ylab=None, xlab=None, title=None):
    if x:
        x = list(x)
        plt.xticks(x)
        plt.bar(x, y)
    else:
        x = list(range(len(y)))
        plt.bar(x, y)
        plt.xticks(x)
    if ylab:
        plt.ylabel(ylab)
    if xlab:
        plt.xlabel(xlab)
    if title:
        plt.title(title)
    plt.show()

def plot_chrom_continuity(vcf_entries, chr_id, x=None, y=None, outfile = None,
                          outfolder = None, returned=False, show=True, label=True, step=1, nb_subplots = None,
                          subplot_init = False, subplot_id = None, title = None, plot_init = False):
    chr_name = list(vcf_entries.keys())[chr_id]
    if label:
        label = chr_name
    if not title:
        title = "Genotyped pos in chr "+str(chr_id+1)+":'"+chr_name+"'"
    chr_entries = vcf_entries[chr_name]
    genotyped_pos = vcf_utils.genotyping_continuity_plot(chr_entries, step=step)
    if returned:
        # if we do not want to plot while executing
        # useful for storing the x,y coords in a variable for ex.
        return genotyped_pos
    else:
        # to plot on the fly
        plot(x=genotyped_pos[0], y=genotyped_pos[1], ylab = "genotyped pos.",
             xlab = "pos. in ref.",
             title = title,
             outfile = outfile, outfolder = outfolder, show=show, label=label,
             nb_subplots = nb_subplots, subplot_init = subplot_init, subplot_id = subplot_id, plot_init = plot_init)

def plot_whole_karyotype(recent_variants, mem_clean = False, step = 1, show = True, min_chr_id = 0,
                         max_chr_id = None, stacked = False, title = None, outfile = None):
    coords = []
    if max_chr_id :
        nb_iter = max_chr_id
    else:
        nb_iter =  len(recent_variants) -1
    if show :
        iter_start = min_chr_id + 1
        if step == "auto" :
            step = round(len(recent_variants[list(recent_variants.keys())[min_chr_id]]) / 1000)
        if stacked:
            nb_subplots = nb_iter - min_chr_id
            subplot_init = True
        else:
            nb_subplots = None
            subplot_init = False
        vcf_utils.customgraphics.plot_chrom_continuity(recent_variants, chr_id = min_chr_id, show = False, returned = False, step = step,
                                                       nb_subplots = nb_subplots, subplot_init = subplot_init, subplot_id = min_chr_id, plot_init = True)
    else :
        iter_start = 0
    for chr in range(iter_start, nb_iter):
        if show == False:
            x, y = vcf_utils.customgraphics.plot_chrom_continuity(recent_variants, chr_id = chr, show = False, returned = True, step = step)
            coords.append([x, y])
            if mem_clean:
                start = time.time()
                del x
                del y
                gc.collect()
                end = time.time()
                print("Cleaned mem. in", str(datetime.timedelta(seconds=end - start)))
        else:
            # if show is enable, use a step
            if step == "auto":
                step = round(len(recent_variants[list(recent_variants.keys())[chr]]) / 1000)
            vcf_utils.customgraphics.plot_chrom_continuity(recent_variants, chr_id = chr, show = False, returned = False, step = step, subplot_id = chr)
        # last case
    if show == True:
        vcf_utils.customgraphics.plot_chrom_continuity(recent_variants, chr_id = nb_iter, show = True, returned = False, step = step, subplot_id = nb_iter,
                                                       title = title,
                                                       outfile = outfile, plot_init = False)
    # maybe add a clean of recent_variants in extreme cases, before building the plots
    if show == False:
        return coords

def plot_chrom_coverage(vcf_entries, chr_id):
    chr_name = list(vcf_entries.keys())[chr_id]
    chr_entries = vcf_entries[chr_name]
    coverage = vcf_utils.compute_coverage(chr_entries)
    barplot(coverage[0], coverage[1], ylab = "coverage (X)",
            xlab = "pos. in ref.",
            title = "Coverage for chr "+str(chr_id+1)+":'"+chr_name+"'")
