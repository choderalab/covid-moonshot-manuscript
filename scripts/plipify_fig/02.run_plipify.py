"""
Assumes the provided glob points to a set of PDB structures with a single ligand / binding site.

Written in the style of asapdiscovery repo scripts `https://github.com/choderalab/covid-moonshot-ml`
"""
import argparse, os, sys
from pathlib import Path
from glob import glob
from tqdm import tqdm
import shutil

from plipify.fingerprints import InteractionFingerprint
from plipify.core import Structure

import pandas as pd

################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    ## Input arguments
    parser.add_argument(
        "-i",
        "--input_glob",
        required=True,
        help="Input glob that will return a list of PDB structures.",
    )

    parser.add_argument("-o", "--output_dir", required=True, help="Path to output directory")

    return parser.parse_args()

def main():
    args = get_args()

    ## Parse symlinks in output_dir
    args.output_dir = os.path.realpath(args.output_dir)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    TARGET = "Mpro"

    HERE = Path(os.getcwd())
    DATA = HERE / "data" / TARGET
    OUT = Path(args.output_dir)
    DATA.mkdir(exist_ok=True, parents=True)
    OUT.mkdir(exist_ok=True, parents=True)

    ## The other structures don't have the same length sequences so bad things happen
    # pdbs = list((DATA / "aligned").glob("Mpro-x*/*_bound_chain*.pdb")) + list((DATA / "aligned").glob("Mpro-z*/*_bound_chain*.pdb"))
    pdbs = list((DATA / "aligned").glob("**/*_bound_chain*.pdb"))
    # print(pdbs)
    ## for debugging
    # pdbs = pdbs[0:200]
    # print(structure_labels)

    ## make new filenames and load structure objects
    structures = []
    for path in tqdm(pdbs, total=len(pdbs)):
        new_fn = OUT / Path(f"{path.parent.name}.pdb")
        if not os.path.exists(new_fn):
            shutil.copy(path, new_fn)
        structure = Structure.from_pdbfile(str(new_fn), ligand_name="LIG")

        ## Skip structures with multiple binding sites
        if len(structure.binding_sites) != 1:
            print(
                f"{path.relative_to(HERE)} contains {len(structure.binding_sites)} binding sites and we want exactly one.")
            continue
        structures.append(structure)
    print(f"Loaded {len(structures)} structures.")

    ## A bunch of filtering
    ## TODO: Use mode instead of median? or only take sequences that have exactly the median?
    ## Remove entries where the difference sequence length - median sequence length is greater than one standard

    ## Collect Sequence and gap penalties
    lengths = pd.DataFrame([((s.identifier), len(s.sequence()), s.sequence()) for s in structures],
                           columns=["identifier", "length", "sequence"])
    lengths["gapcount"] = lengths.sequence.str.count('-')

    print('Sequence length median and std: ', lengths.length.median(), lengths.length.std())
    print(f"{sum(lengths.gapcount > 0)} sequences with gaps to be removed.")


    ## Perform the actual filtering
    lengths = lengths[lengths.gapcount == 0]
    # lengths = lengths[(lengths.length - lengths.length.median()).abs() < lengths.length.std()]
    lengths = lengths[lengths.length == lengths.length.median()]
    filtered_structures = [s for s in structures if s.identifier in set(lengths.identifier.tolist())]
    print(len(pdbs), "->", len(structures), "->", len(filtered_structures), "=", len(pdbs) - len(filtered_structures),
          "structures filtered out")

    structures = filtered_structures

    structure_name_type_dict = {'Mpro-P': 0, 'Mpro-x': 0, 'Mpro-z': 0, 'other': 0}
    print(structure_name_type_dict.keys())
    for s in structures:
        name_code = s.identifier[0:6]
        if name_code in structure_name_type_dict.keys():
            structure_name_type_dict[name_code] += 1
        else:
            structure_name_type_dict['other'] += 1
    print(structure_name_type_dict)

    ## LEGACY CODE, LEAVING FOR NOW
    # # exclude 35 Mpro-P structures
    # filtered_structures = [s for s in structures if not (s.identifier.startswith('Mpro-J'))]
    # print('remaining: ', len(filtered_structures), 'from:', len(structures))
    # structures = filtered_structures

    # Review
    fp = InteractionFingerprint().calculate_fingerprint(
        structures,
        labeled=True,
        as_dataframe=True,
        remove_non_interacting_residues=True,
        remove_empty_interaction_types=True,
        ensure_same_sequence=False
    )
    out_csv = os.path.join(args.output_dir, "plipify_results.csv")
    fp.to_csv(out_csv)

if __name__ == "__main__":
    main()