from collections import Counter, OrderedDict

import numpy as np


def rxn_order_corr(rxn_correlation):
    # calculate the absolute matrix so we can find maximum correlation values
    abs_rxn_corr = np.absolute(rxn_correlation)

    max_inds = []
    correlation_inds = {}
    for i in range(len(rxn_correlation)):
        highly_correlated_met = np.argmax(np.absolute(rxn_correlation[i]))
        max_inds.append(np.argmax(highly_correlated_met))

        if highly_correlated_met in correlation_inds:
            correlation_inds[highly_correlated_met].append(i)
        else:
            correlation_inds[highly_correlated_met] = [i]

    counts = Counter(max_inds)
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    order = [sorted_counts[i][0] for i in range(len(sorted_counts))]

    finished = OrderedDict()
    for ele in order:
        if ele in finished.keys():
            continue
        for i in range(1, len(rxn_correlation) + 1):
            largest_cor = np.partition(rxn_correlation[:, ele].flatten(), -i)[-i]
            ind_largest_cor = np.where(
                largest_cor == np.absolute(rxn_correlation[:, ele])
            )
            if ind_largest_cor in finished:
                finished[ele] = ind_largest_cor
        # if still can't fid suitable warm start ignore it
        finished[ele] = "NA"

        # Now add the elements related to ele
        if ele in correlation_inds:
            for item in correlation_inds[ele]:
                finished[item] = ele
