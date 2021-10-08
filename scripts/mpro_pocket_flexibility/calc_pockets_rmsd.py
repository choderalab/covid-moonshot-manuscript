import argparse

import matplotlib.pyplot as plt
import MDAnalysis as mda
import pandas as pd
import seaborn as sns
from MDAnalysis.analysis.rms import rmsd
from numpy.core.fromnumeric import amin

parser = argparse.ArgumentParser(description="Calculate MPro pocket flexibility")
parser.add_argument(
    "-path",
    dest="path",
    type=str,
    help="the path to the fragment directories e.g. path/to/fragments/aligned/",
)
parser.add_argument(
    "-metadata_file",
    dest="metadata_file",
    type=str,
    help="the metadata.csv file with relevant fragment information",
)

args = parser.parse_args()


def calc_rmsd(
    fragments: list,
    pockets: dict,
    path: str,
    superposition: bool = False,
    backbone: bool = False,
) -> dict:
    results = {fragment: {} for fragment in fragments}

    for fragment in fragments:
        for pocket in pockets:

            # Create MDAnalysis Universe object for current fragment PDB
            u = mda.Universe(path + f"{fragment}/" + f"{fragment}_apo-desolv.pdb")

            # make selections, don't count alternative locations (?)
            if backbone:
                pocket_sel = (
                    "backbone and resid " + f"{pockets[pocket]} and not altloc B"
                )
                ref_pocket_sel = "backbone and resid " + f"{pockets[pocket]}"

            else:
                pocket_sel = (
                    "not backbone and not name H* and resid "
                    + f"{pockets[pocket]} and not altloc B"
                )
                ref_pocket_sel = (
                    "not backbone and not name H*  and resid " + f"{pockets[pocket]}"
                )

            # Calculate RMSD of current pocket selection against reference pocket
            pocket_rmsd = rmsd(
                u.select_atoms(pocket_sel).positions,
                ref_u.select_atoms(ref_pocket_sel).positions,
                superposition=superposition,
            )

            # Populate dictionary with results for each fragment and each pocket
            results[fragment][pocket] = pocket_rmsd

    return results


def sort_results(data: dict, pocket_string: str) -> dict:

    # dict_to_sort = data.pop("series")
    dict_to_sort = data

    sorted_keys = sorted(dict_to_sort, key=lambda x: (dict_to_sort[x][pocket_string]))

    return sorted_keys


def plot_pocket_hist(df: pd.DataFrame, save_name: str, legend=False, despine_y=False) -> None:

    pymol_colours = [
        (0.99, 0.82, 0.65), # amino - wheat
        (0.65, 0.9, 0.65), # ugi - pale green
        (1.0, 0.5, 1.0), # quin - violet
        (0.96, 0.41, 0.23), # benzo - cb_orange
    ]
    
    sns.set_palette(sns.color_palette(pymol_colours))

    fig, ax = plt.subplots(figsize = (1.5,1.5))

    g = sns.kdeplot(
        data=df,
        x="RMSD",
        hue="Series",
        fill=True,
        alpha=0.75,
        linewidth=1,
        common_norm=False,
        legend=False,
        ax=ax,
    )

    # control figure aesthetics
    g.set(yticklabels=[])
    g.set(yticks=[])
    if not despine_y:
        g.set(ylabel="Number of structures")
    g.set(xlabel=r"Pocket RMSD ($\AA$)")
    g.set(xlim=[0, 4])
    g.set(title=df.Pocket.tolist()[0])

    # Remove top and right box lines
    sns.despine()
    if despine_y:
        sns.despine(left=True)
        g.set(ylabel="")

    # This is a hack as the legend isn't directly accessible from g
    if legend:
        plt.legend(
            labels = ['Benzotriazoles','Quinolones', 'Ugis', 'Aminopyridines'],
            loc = 'right',
            fontsize=4.5,
            frameon=False,
            )

    plt.tight_layout()

    print(f"Saving plot as {save_name}.png")
    plt.savefig(f"{save_name}_hist.png", dpi=300)

    plt.clf()


def get_fragment_list(metadata_file: str, series: str) -> list:

    series_dict = {
        "amino": "Aminopyridine-like",
        "ugi": "Ugi",
        "quin": "Quinolone",
        "benzo": "Benzotriazole",
    }

    # Read metadata and create dataframe
    df = pd.read_csv(metadata_file)

    # Create dataframe for specified series
    df_series = df[df["site_name"] == series_dict[series]]

    # Get the crystal names
    df_series_xtal_names = df_series["crystal_name"].to_list()

    df_series_xtal_names

    return df_series_xtal_names


def get_pocket_rmsd(pocket: str, result: dict, series: str):

    full_series_names = {
        "amino": "Aminopyridines",
        "ugi": "Ugis",
        "quin": "Quinolones",
        "benzo": "Benzotriazoles"

    }

    full_pocket_names = {
        "P1": "P1",
        "P1_prime": r"P1$^\prime$",
        "P2": "P2",
        "P3_4_5": "P3-5"

    }

    rmsd_results = {"RMSD": [], "Pocket": full_pocket_names[pocket], "Series": full_series_names[series]}

    for fragment in result:
        rmsd_results["RMSD"].append(result[fragment][pocket])

    df = pd.DataFrame(rmsd_results)

    return df


if __name__ == "__main__":

    path = args.path
    metadata_file = args.metadata_file

    # methoxy benzo ref = Mpro-P0157_0A
    ref_code = "Mpro-P0157_0A"
    ref_prot = path + f"{ref_code}/{ref_code}_apo-desolv.pdb"
    ref_u = mda.Universe(ref_prot)

    # Define pocket residues
    pockets = {
        "P1": "142 141 140 172 163 143 144",
        "P1_prime": "25 26 27",
        "P2": "41 49 54",
        "P3_4_5": "189 190 191 192 168 167 166 165",
    }

    # TODO put this in a for loop
    amino_fragments = get_fragment_list(metadata_file, "amino")
    ugi_fragments = get_fragment_list(metadata_file, "ugi")
    quin_fragments = get_fragment_list(metadata_file, "quin")
    benzo_fragments = get_fragment_list(metadata_file, "benzo")

    # Calculate the RMSD for each series
    amino_result = calc_rmsd(amino_fragments, pockets, path)
    ugi_result = calc_rmsd(ugi_fragments, pockets, path)
    quin_result = calc_rmsd(quin_fragments, pockets, path)
    benzo_result = calc_rmsd(benzo_fragments, pockets, path)

    # Create large dataframe for histogram plots
    # these are RMSDs of each pocket with different
    # histograms for each lead series

    d = {
        "amino": amino_result,
        "ugi": ugi_result,
        "quin": quin_result,
        "benzo": benzo_result,
    }

    for pocket in ["P1", "P1_prime", "P2", "P3_4_5"]:
        legend = False
        df_store = []
        for series, result in d.items():

            df_store.append(get_pocket_rmsd(pocket, result, series))

        df = pd.concat(df_store)

        # Plotting
        if pocket == "P1_prime": # add legend to one plot
            legend = True
        plot_pocket_hist(df, pocket, legend=legend, despine_y=True)

    # Get a sorted list (lowest to highest) of RMSDs
    amino_p1_displacement = sort_results(amino_result, "P1")
    amino_p1_prime_displacement = sort_results(amino_result, "P1_prime")
    amino_p2_displacement = sort_results(amino_result, "P2")
    amino_p3_4_5_displacement = sort_results(amino_result, "P3_4_5")

    ugi_p1_displacement = sort_results(ugi_result, "P1")
    ugi_p1_prime_displacement = sort_results(ugi_result, "P1_prime")
    ugi_p2_displacement = sort_results(ugi_result, "P2")
    ugi_p3_4_5_displacement = sort_results(ugi_result, "P3_4_5")

    quin_p1_displacement = sort_results(quin_result, "P1")
    quin_p1_prime_displacement = sort_results(quin_result, "P1_prime")
    quin_p2_displacement = sort_results(quin_result, "P2")
    quin_p3_4_5_displacement = sort_results(quin_result, "P3_4_5")

    benzo_p1_displacement = sort_results(benzo_result, "P1")
    benzo_p1_prime_displacement = sort_results(benzo_result, "P1_prime")
    benzo_p2_displacement = sort_results(benzo_result, "P2")
    benzo_p3_4_5_displacement = sort_results(benzo_result, "P3_4_5")

    amino_rmsd = [
        amino_result[amino_p1_displacement[-1]]["P1"],
        amino_result[amino_p1_prime_displacement[-1]]["P1_prime"],
        amino_result[amino_p2_displacement[-1]]["P2"],
        amino_result[amino_p3_4_5_displacement[-1]]["P3_4_5"],
    ]

    ugi_rmsd = [
        ugi_result[ugi_p1_displacement[-1]]["P1"],
        ugi_result[ugi_p1_prime_displacement[-1]]["P1_prime"],
        ugi_result[ugi_p2_displacement[-1]]["P2"],
        ugi_result[ugi_p3_4_5_displacement[-1]]["P3_4_5"],
    ]

    quin_rmsd = [
        quin_result[quin_p1_displacement[-1]]["P1"],
        quin_result[quin_p1_prime_displacement[-1]]["P1_prime"],
        quin_result[quin_p2_displacement[-1]]["P2"],
        quin_result[quin_p3_4_5_displacement[-1]]["P3_4_5"],
    ]

    benzo_rmsd = [
        benzo_result[benzo_p1_displacement[-1]]["P1"],
        benzo_result[benzo_p1_prime_displacement[-1]]["P1_prime"],
        benzo_result[benzo_p2_displacement[-1]]["P2"],
        benzo_result[benzo_p3_4_5_displacement[-1]]["P3_4_5"],
    ]

    max_df_rmsd = pd.DataFrame(
        {
            "amino": amino_rmsd,
            "ugi": ugi_rmsd,
            "quin": quin_rmsd,
            "benzo": benzo_rmsd,
        }
    )

    max_df = pd.DataFrame(
        {
            "amino": [
                amino_p1_displacement[-1],
                amino_p1_prime_displacement[-1],
                amino_p2_displacement[-1],
                amino_p3_4_5_displacement[-1],
            ],
            "ugi": [
                ugi_p1_displacement[-1],
                ugi_p1_prime_displacement[-1],
                ugi_p2_displacement[-1],
                ugi_p3_4_5_displacement[-1],
            ],
            "quin": [
                quin_p1_displacement[-1],
                quin_p1_prime_displacement[-1],
                quin_p2_displacement[-1],
                quin_p3_4_5_displacement[-1],
            ],
            "benzo": [
                benzo_p1_displacement[-1],
                benzo_p1_prime_displacement[-1],
                benzo_p2_displacement[-1],
                benzo_p3_4_5_displacement[-1],
            ],
        }
    )

    # Sort out the row names
    max_df = max_df.rename(index={0: "P1", 1: "P1_prime", 2: "P2", 3: "P3_4_5"})
    max_df_rmsd = max_df_rmsd.rename(
        index={0: "P1", 1: "P1_prime", 2: "P2", 3: "P3_4_5"}
    )

    print(max_df)
    print(max_df_rmsd)

    # Save CSV files, one for fragment names and one for RMSD values
    csv_name = "mpro_pockets_max_flex_fragments.csv"
    print(f"writing fragment results to {csv_name}")
    max_df.to_csv(csv_name)

    csv_name2 = "mpro_pockets_max_flex_rmsds.csv"
    print(f"writing rmsd values to {csv_name2}")
    max_df_rmsd.to_csv(csv_name2)
