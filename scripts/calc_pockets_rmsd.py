import argparse
import MDAnalysis as mda
import pandas as pd
import os

from MDAnalysis.analysis.rms import rmsd

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


def calc_rmsd(fragments, pockets, path, superposition=False, backbone=False):
    df = {fragment: {} for fragment in fragments}

    for fragment in fragments:
        for pocket in pockets:

            u = mda.Universe(path + f"{fragment}/" + f"{fragment}_apo-desolv.pdb")

            # make selections, don't count alternative locations (?)
            if backbone:
                pocket_sel = (
                    "backbone and resid " + f"{pockets[pocket]} and not altloc B"
                )
                ref_pocket_sel = "backbone and resid " + f"{pockets[pocket]}"

            else:
                pocket_sel = (
                    "not name H* and resid " + f"{pockets[pocket]} and not altloc B"
                )
                ref_pocket_sel = "not name H*  and resid " + f"{pockets[pocket]}"

            # Calculate RMSD of current pocket selection against reference pocket
            pocket_rmsd = rmsd(
                u.select_atoms(pocket_sel).positions,
                ref_u.select_atoms(ref_pocket_sel).positions,
                superposition=superposition,
            )

            df[fragment][pocket] = pocket_rmsd

    return df


def sort_results(data, pocket_string):

    sorted_keys = sorted(data, key=lambda x: (data[x][pocket_string]))

    return sorted_keys


def get_fragment_list(metadata_file, series):

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

    return df_series_xtal_names


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

    amino_result = calc_rmsd(amino_fragments, pockets, path)
    ugi_result = calc_rmsd(ugi_fragments, pockets, path)
    quin_result = calc_rmsd(quin_fragments, pockets, path)
    benzo_result = calc_rmsd(benzo_fragments, pockets, path)

    # get a sorted list (highest to lowest) of RMSDs
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

    print(amino_result[amino_p1_displacement[0]])
    print(amino_result[amino_p1_displacement[0]]["P1"])

    amino_rmsd = [
        amino_result[amino_p1_displacement[0]]["P1"],
        amino_result[amino_p1_prime_displacement[0]]["P1_prime"],
        amino_result[amino_p2_displacement[0]]["P2"],
        amino_result[amino_p3_4_5_displacement[0]]["P3_4_5"],
    ]

    ugi_rmsd = [
        ugi_result[ugi_p1_displacement[0]]["P1"],
        ugi_result[ugi_p1_prime_displacement[0]]["P1_prime"],
        ugi_result[ugi_p2_displacement[0]]["P2"],
        ugi_result[ugi_p3_4_5_displacement[0]]["P3_4_5"],
    ]

    quin_rmsd = [
        quin_result[quin_p1_displacement[0]]["P1"],
        quin_result[quin_p1_prime_displacement[0]]["P1_prime"],
        quin_result[quin_p2_displacement[0]]["P2"],
        quin_result[quin_p3_4_5_displacement[0]]["P3_4_5"],
    ]

    benzo_rmsd = [
        benzo_result[benzo_p1_displacement[0]]["P1"],
        benzo_result[benzo_p1_prime_displacement[0]]["P1_prime"],
        benzo_result[benzo_p2_displacement[0]]["P2"],
        benzo_result[benzo_p3_4_5_displacement[0]]["P3_4_5"],
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
                amino_p1_displacement[0],
                amino_p1_prime_displacement[0],
                amino_p2_displacement[0],
                amino_p3_4_5_displacement[0],
            ],
            "ugi": [
                ugi_p1_displacement[0],
                ugi_p1_prime_displacement[0],
                ugi_p2_displacement[0],
                ugi_p3_4_5_displacement[0],
            ],
            "quin": [
                quin_p1_displacement[0],
                quin_p1_prime_displacement[0],
                quin_p2_displacement[0],
                quin_p3_4_5_displacement[0],
            ],
            "benzo": [
                benzo_p1_displacement[0],
                benzo_p1_prime_displacement[0],
                benzo_p2_displacement[0],
                benzo_p3_4_5_displacement[0],
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
    max_df.to_csv(csv_name2)
