import numpy as np
from collections import OrderedDict, Counter


def rxn_order_corr(rxn_correlation):
    # calculate the absolute matrix so we can find maximum correlation values
    abs_rxn_corr = np.absolute(rxn_correlation)

    max_inds = []
    for i in range(len(rxn_correlation)):
        max_inds.append(np.argmax(np.absolute(rxn_correlation)))

    counts = Counter(max_inds)
    sorted_counts = sorted(counts.items(), key=lambda x: x[1], reverse=True)
    order = [sorted_counts[i][0] for i in range(len(sorted_counts))]

    finished = OrderedDict()
    for ele in order:
        if ele in finished.keys():
            continue

        # Find max correlated for ele is present in finished
        finished[ele] = "NA"
        # Find all rxn indices with ele as highly correlated
        matched_inds = [ind for ind, val in enumerate(max_inds) if val == ele]

        for ind in matched_inds:
            finished[ind] = ele

    warm_start = OrderedDict()
    finished = []
    for i in range(len(rxn_correlation)):
        if i in finished:
            continue
        temp = rxn_correlation[:, i]
        for j in range(1, len(rxn_correlation) + 1):
            best_correlation = np.where(
                np.absolute(temp) == np.partition(np.absolute(temp).flatten(), -j)[-j]
            )[0][0]
            if best_correlation not in finished:
                continue
            warm_start[i] = best_correlation
        finished.append(i)
