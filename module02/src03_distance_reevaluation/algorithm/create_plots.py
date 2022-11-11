import matplotlib.pyplot as plt
import pandas as pd
import numpy as np


def create_histogram(data, filename, output_directory):
    # Create histogram plots with matplotlib
    #
    # input data: pd: DataFrame, bins: list(int), filename: str, output_directory: str
    # no return

    # normal distance histograms
    dist_data = data[["eucl_dist_tplk", "topo_dist_tplk"]].copy()
    bins = [x * 10 for x in range(9)]
    bin_centers = np.diff(bins) * .5 + bins[:-1]
    dist_data[dist_data > 70] = 70
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    known_xy = []
    for col in dist_data.columns:
        freq, _, patches = plt.hist(dist_data[[col]], bins=bins, alpha=.3,
                                    label=f"{col} (n={len(dist_data[~pd.isna(dist_data[col])].index)})")
        i = 0
        for fr, x, patch in zip(freq, bin_centers, patches):
            height = int(freq[i])
            y = height
            if x in [val for val, _ in known_xy]:
                known_y = known_xy[[val for val, _ in known_xy].index(x)][1]
                rad = int(round(max(int(f) for f in freq) / 40) + 1)
                if y in range(known_y-rad, known_y+rad+1) and y != known_y:
                    y = known_y - rad if y <= known_y else known_y + rad
            plt.annotate(str(height) if height > 0 else '', xy=(x, y), xytext=(0, .2),
                         textcoords="offset points", ha="center", va="bottom")
            known_xy.append((x, y))
            i += 1
    plt.xlabel(f"distance " + r"[$\AA$]")
    plt.ylabel("frequency")
    x_labels = [str(x) for x in bins]
    x_labels[-1] = str(float("inf"))
    plt.xticks(bins, x_labels)
    plt.legend()
    plt.title("Frequency of interaction site distances in structures")
    plt.savefig(f"{output_directory}{filename}_hist.png")

    # Clear figure
    plt.clf()

    # difference distance histograms
    diff_data = pd.DataFrame()
    diff_data["self_eucl - tplk_topo"] = abs(data["eucl_dist"] - data["topo_dist_tplk"])
    diff_data["tplk_eucl - tplk_topo"] = abs(data["eucl_dist_tplk"] - data["topo_dist_tplk"])
    diff_data[diff_data > 11] = 11
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    known_xy = []
    if not pd.isna(diff_data.max().max()):
        bins = list(range(12))
        bin_centers = np.diff(bins) * .5 + bins[:-1]
        for col in diff_data.columns:
            freq, _, patches = plt.hist(diff_data[col], bins=bins, alpha=.3,
                                        label=f"{col} (n={len(diff_data[~pd.isna(diff_data[col])].index)})")
            i = 0
            for fr, x, patch in zip(freq, bin_centers, patches):
                height = int(freq[i])
                y = height
                if x in [val for val, _ in known_xy]:
                    known_y = known_xy[[val for val, _ in known_xy].index(x)][1]
                    rad = int(round(max(int(f) for f in freq) / 40) + 1)
                    if y in range(known_y-rad, known_y+rad+1) and y != known_y:
                        y = known_y - rad if y <= known_y else known_y + rad
                plt.annotate(str(height) if height > 0 else '', xy=(x, y), xytext=(0, .2),
                             textcoords="offset points", ha="center", va="bottom")
                known_xy.append((x, y))
                i += 1
        plt.xlabel(f"absolute distance difference " + r"[$\AA$]")
        plt.ylabel("frequency")
        x_labels = [str(x) for x in bins]
        x_labels[-1] = "inf"
        plt.xticks(bins, x_labels)
        plt.legend()
        plt.title("Frequency of absolute interaction site distance differences")
        plt.savefig(f"{output_directory}{filename}_hist_diff.png")

    # Clear figure
    plt.clf()

    # pdb method histogram
    method_data = data[["pdb_method", "best_res_pdb_method"]].copy().astype(str)
    label_dict = {x: i for i, x in enumerate(sorted(pd.unique(method_data.values.ravel('K'))))}
    bins = list(range(len(label_dict) + 1))
    bin_centers = np.diff(bins) * .5 + bins[:-1]
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    known_xy = []
    for col in method_data.columns:
        freq, _, patches = plt.hist(method_data[col].map(label_dict), bins=bins, alpha=.3,
                                    label=f"{'final' if col.startswith('pdb') else 'best struct search res'}")
        i = 0
        for fr, x, patch in zip(freq, bin_centers, patches):
            height = int(freq[i])
            y = height
            if x in [val for val, _ in known_xy]:
                known_y = known_xy[[val for val, _ in known_xy].index(x)][1]
                rad = int(round(max(int(f) for f in freq) / 40) + 1)
                if y in range(known_y-rad, known_y+rad+1) and y != known_y:
                    y = known_y - rad if y <= known_y else known_y + rad
            plt.annotate(str(height) if height > 0 else '', xy=(x, y), xytext=(0, .2),
                         textcoords="offset points", ha="center", va="bottom")
            known_xy.append((x, y))
            i += 1
    plt.xlabel(f"pdb experimental method")
    plt.ylabel("frequency")
    x_labels = [x.replace(' ', '\n') if x != "nan" else "unknown"
                for x in sorted(pd.unique(method_data.values.ravel('K')))]
    x_labels = [x if x != '-' else "not found" for x in x_labels]
    plt.xticks(bin_centers, x_labels, rotation=45)
    plt.legend()
    plt.title("Frequency of experimental methods in structures")
    plt.savefig(f"{output_directory}{filename}_pdb_exp.png")

    # Clear figure
    plt.clf()

    # pdb resolution histogram
    res_data = data[["pdb_resolution", "best_res_pdb_resolution"]].copy()
    res_data[res_data == "ALPHAFOLD"] = -1
    res_data[res_data == '-'] = -2
    res_data = res_data.astype(float)
    res_data[res_data > 13] = 14
    bins = [x - 2 for x in range(18)]
    bin_centers = np.diff(bins) * .5 + bins[:-1]
    plt.figure(figsize=(6.5, 6), constrained_layout=True)
    known_xy = []
    for col in res_data.columns:
        freq, _, patches = plt.hist(res_data[col], bins=bins, alpha=.3,
                                    label=f"{'final' if col.startswith('pdb') else 'best struct search res'}")
        i = 0
        for fr, x, patch in zip(freq, bin_centers, patches):
            height = int(freq[i])
            y = height
            if x in [val for val, _ in known_xy]:
                known_y = known_xy[[val for val, _ in known_xy].index(x)][1]
                rad = int(round(max(int(f) for f in freq) / 40) + 1)
                if y in range(known_y-rad, known_y+rad+1) and y != known_y:
                    y = known_y - rad if y <= known_y else known_y + rad
            plt.annotate(str(height) if height > 0 else '', xy=(x, y), xytext=(0, .2),
                         textcoords="offset points", ha="center", va="bottom")
            known_xy.append((x, y))
            i += 1
    plt.xlabel(f"pdb resolution")
    plt.ylabel("frequency")
    x_labels = [f"[{x},{x+1})" if x != bins[-2] else f"[{x},{x+1}]" for x in bins[:-1]]
    x_labels[1] = "ALPHA-\nFOLD"
    x_labels[0] = "not\nfound"
    plt.xticks(bin_centers, x_labels, rotation=45, fontsize=8)
    plt.legend()
    plt.title("Frequency of resolutions in structures")
    plt.savefig(f"{output_directory}{filename}_pdb_res.png")
