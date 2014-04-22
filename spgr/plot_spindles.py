from matplotlib.pyplot import hist, plot, imshow, subplots


def hist_values(all_subj, all_spindles,
                get_value, x_lim, nbin):

    f, subp = subplots(3, 2, figsize=(18, 12))

    for ax, subj, spindles in zip(subp.flatten(), all_subj, all_spindles):
        value = [get_value(x) for x in spindles]
        ax.hist(value, range=x_lim,
                bins=(x_lim[1] - x_lim[0]) * nbin,
                align='left')
        ax.set_title(subj)
        ax.set_xlim(x_lim)
