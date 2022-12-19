"""
Assumes the provided glob points to a set of PDB structures with a single ligand / binding site.

Written in the style of asapdiscovery repo scripts `https://github.com/choderalab/covid-moonshot-ml`
"""
import argparse, os, sys
from plipify.fingerprints import InteractionFingerprint

################################################################################
def get_args():
    parser = argparse.ArgumentParser(description="")

    ## Input arguments
    parser.add_argument(
        "-i",
        "--input_glob",
        help="Input glob that will return a list of PDB structures.",
    )
    parser.add_argument(
        "-n",
        "--n_structures",
        default=None,
        help="How many of the provided structures to run plipify on. "
             "Defaults to None, which means all structures. Useful for debugging."
    )


    return parser.parse_args()


def main():
    args = get_args()

    ## Parse symlinks in output_dir
    args.output_dir = os.path.realpath(args.output_dir)

    # Review
    fp = InteractionFingerprint().calculate_fingerprint(
        structures,
        labeled=True,
        as_dataframe=True,
        remove_non_interacting_residues=True,
        remove_empty_interaction_types=True,
        ensure_same_sequence=False
    )


if __name__ == "__main__":
    main()