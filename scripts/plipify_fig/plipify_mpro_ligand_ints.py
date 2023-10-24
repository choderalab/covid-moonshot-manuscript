import os
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile

import pandas as pd
import plipify.core
import plipify.fingerprints
import plipify.visualization
import requests
from Bio.PDB import *
from plipify.core import Structure
from plipify.fingerprints import InteractionFingerprint
from plipify.visualization import PymolVisualizer, fingerprint_writepdb
from tqdm.auto import tqdm

TARGET = "Mpro"

HERE = Path(os.getcwd())
DATA = HERE / "data" / TARGET
OUT = HERE / "output" / TARGET
DATA.mkdir(exist_ok=True, parents=True)
OUT.mkdir(exist_ok=True, parents=True)

r = requests.get(f"https://fragalysis.diamond.ac.uk/api/targets/?format=json&title={TARGET}")
r.raise_for_status()
target = r.json()["results"][0]

if target["zip_archive"] is None:
    raise ValueError(f"Target {TARGET} is not downloadable")

# If already downloaded, they should be here:
pdbs = list((DATA / "aligned").glob("**/*_bound.pdb"))

if not pdbs:
    archive =  requests.get(target["zip_archive"], stream=True)
    with BytesIO(archive.raw.data) as b, ZipFile(b) as z:
        z.extractall(DATA)
    # Reassign now
    pdbs = list((DATA / "aligned").glob("**/*_bound.pdb"))
    assert pdbs, "Couldn't find downloaded PDB structures!"

for nr, filepath in enumerate(pdbs):
    pdb_id = str(pdbs[nr]).split('/')[-1][:-4]
    chain_id = str(pdbs[nr]).split('/')[-1].split('_')[1][-1]
    new_filename=str(pdbs[nr])[:-4]+'_chain'+str(chain_id)+'.pdb'

    ## Read the PDB file and extract the chain from structure[0]
    model = PDBParser(PERMISSIVE=1,QUIET=1).get_structure(pdb_id, filepath)[0]
    ### Save new file
    io = PDBIO()
    io.set_structure(model[chain_id])
    io.save(new_filename)

# Update the pdb list to pdb chain list
# pdbs = list((DATA / "aligned").glob("**/*_bound_chain*.pdb"))

## The other structures don't have the same length sequences so bad things happen
pdbs = list((DATA / "aligned").glob("Mpro-x*/*_bound_chain*.pdb")) + list((DATA / "aligned").glob("Mpro-z*/*_bound_chain*.pdb"))

# Use plipify
structures = []
for path in tqdm(pdbs):
    structure = Structure.from_pdbfile(str(path), ligand_name="LIG")
    if len(structure.binding_sites) != 1:
        print(f"{path.relative_to(HERE)} contains {len(structure.binding_sites)} binding sites and we want exactly one.")
        continue
    structures.append(structure)

lengths = pd.DataFrame([((s.identifier), len(s.sequence())) for s in structures], columns=["identifier", "length"])
# Remove entries where the difference sequence length - median sequence length is greater than one standard deviation
print('Sequence length median and std: ',lengths.length.median(),lengths.length.std() )

lengths = lengths[(lengths.length - lengths.length.median()).abs() < lengths.length.std()]
filtered_structures = [s for s in structures if s.identifier in set(lengths.identifier.tolist())]
print(len(pdbs), "->", len(structures), "->", len(filtered_structures), "=", len(pdbs) - len(filtered_structures), "structures filtered out")

structures = filtered_structures

structure_name_type_dict={'Mpro-P':0, 'Mpro-x':0, 'Mpro-z':0, 'other':0}
print(structure_name_type_dict.keys())
for s in structures:
    name_code = s.identifier[0:6]
    if name_code in structure_name_type_dict.keys():
        structure_name_type_dict[name_code]+=1
    else:
        structure_name_type_dict['other']+=1
print(structure_name_type_dict)

# exclude 35 Mpro-P structures
filtered_structures = [s for s in structures if not (s.identifier.startswith('Mpro-P'))]
print('remaining: ',len(filtered_structures), 'from:', len(structures))
structures = filtered_structures

# Review
fp = InteractionFingerprint().calculate_fingerprint(
        structures,
        labeled=True, 
        as_dataframe=True, 
        remove_non_interacting_residues=True,
        remove_empty_interaction_types=True,
        ensure_same_sequence=False
    )


if not fp.values.shape[0]:
    raise ValueError("Fingerprint is empty!")

print(fp)

# filter out some of the low counts
fp_focused = fp[fp.sum(axis=1) > 10]

print(fp_focused)

# make images
methoxy = Structure.from_pdbfile(path="./data/Mpro/aligned/Mpro-P0157_0A/Mpro-P0157_0A_bound_chainA.pdb")

output_path = Path("./output/")
pdb_ints = fingerprint_writepdb(fingerprint_df=fp_focused, structure=methoxy, output_path=output_path, ligand=True, summed=True)

spectrum_cols = ["white_green", "white_pink", "white_cyan", "white_red", "white_blue", "white_yellow", "white_orange", "white_purple"]

for i, interaction in enumerate(pdb_ints):

    v = PymolVisualizer(pdb=pdb_ints[interaction]._path)
    v.set_style()
    v.create_image(surface=True, ligand_col="cyan", spectrum_col=spectrum_cols[i], show_ligand=False, cnc_protein=True)
    v.render(name=f"{interaction}_pymol_image", save_path="./output/")
