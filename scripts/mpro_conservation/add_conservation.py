import argparse
import os

import MDAnalysis as mda
import pandas as pd

parser = argparse.ArgumentParser(description="Calculate MPro pocket flexibility")
parser.add_argument(
    "-pdb_file",
    dest="pdb_file",
    type=str,
    help="the input PDB file",
)
parser.add_argument(
    "-alignment_file",
    dest="alignment_file",
    type=str,
    help="an alignment CSV file containing conservation information",
)

args = parser.parse_args()


def assign_conservation(
    pdb_file: str,
    alignment_file: str,
) -> None:
    def _read_csv(alignment_file: str) -> dict:
        df = pd.read_csv(alignment_file)

        # Extract relevant information
        df_bfactors = df[["Res Num", "Identity"]]
        df_bfactors.columns = ["resnum", "id"]

        resid_bfactors = dict(zip(df_bfactors.resnum, df_bfactors.id))

        return resid_bfactors

    u = mda.Universe(pdb_file)

    # make all bfactor values zero
    u.add_TopologyAttr("tempfactors")

    # get a dict of resids and bfactors to assign
    resid_bfactors = _read_csv(alignment_file)

    # assign bfactors to every atom in a given residue
    for resid in resid_bfactors:
        sel = u.select_atoms(f"resid {resid}")
        sel.tempfactors = resid_bfactors[resid]

    # write out the final structure
    output_name = f"{os.path.basename(pdb_file).split('.')[0]}_conservation.pdb"
    with mda.Writer(output_name, n_atoms=u.atoms.n_atoms) as PDB:
        PDB.write(u.atoms)

    print(f"Wrote PDB file {output_name} with conservation in bfactor field")


if __name__ == "__main__":

    pdb_file = args.pdb_file
    alignment_file = args.alignment_file

    assign_conservation(pdb_file, alignment_file)
