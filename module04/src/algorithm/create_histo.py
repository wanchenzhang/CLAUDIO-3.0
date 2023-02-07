import matplotlib.pyplot as plt
import numpy as np


def create_histograms(data, intra_only, filename, cutoff, output_directory):
    # create histograms for inter scores
    #
    # input data: pd.DataFrame, intra_only: bool, filename: str, output_directory: str
    # no return

    if intra_only:
        inter_data = data[data.inter_score > cutoff][["inter_score"]]
    else:
        inter_data = data[(data.inter_score > cutoff) & (data.unip_id_a == data.unip_id_b)][["inter_score"]]
    if not inter_data.empty:
        label = "experimental confidence score for intra crosslinks (" \
                + "$n_{>cutoff}$" + f" = {sum(inter_data.inter_score > cutoff)})"
        bins = np.array(list(range(11))) * .1
        bin_centers = np.diff(bins) * .5 + bins[:-1]
        plt.figure(figsize=(6.5, 6), constrained_layout=True)
        known_xy = []
        freq, _, _ = plt.hist(inter_data, bins=bins, alpha=.3, label=label)
        for fr, x in zip(freq, bin_centers):
            height = int(fr)
            y = height
            if x in [val for val, _ in known_xy]:
                known_y = known_xy[[val for val, _ in known_xy].index(x)][1]
                rad = int(round(max(int(f) for f in freq) / 40) + 1)
                if y in range(known_y - rad, known_y + rad + 1) and y != known_y:
                    y = known_y - rad if y <= known_y else known_y + rad
            plt.annotate(str(height) if height > 0 else '', xy=(x, y), xytext=(0, .2), textcoords="offset points",
                         ha="center", va="bottom")
            known_xy.append((x, y))
        plt.xlabel(f"score")
        plt.ylabel("frequency")
        x_labels = [f"{x:.1f}" for x in bins]
        plt.xticks(bins, x_labels)
        plt.title(f"Inter Interaction Score (Reclassification cutoff = {cutoff})")
        plt.legend()
        plt.savefig(f"{output_directory}{filename}_inter_score.png")
