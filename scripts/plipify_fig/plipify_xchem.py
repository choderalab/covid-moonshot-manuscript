import os

import requests
r = requests.get(f"https://fragalysis.diamond.ac.uk/api/targets/?format=json")
r.raise_for_status()
targets = [t["title"] for t in r.json()["results"] if t["zip_archive"]]
print("Downloadable targets:", *targets)

from pathlib import Path

TARGET = "Mpro"  # other options, see above

#HERE = Path(_dh[-1])
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

from io import BytesIO
from zipfile import ZipFile

# If already downloaded, they should be here:
pdbs = list((DATA / "aligned").glob("**/*_bound.pdb"))

if not pdbs:
    archive =  requests.get(target["zip_archive"], stream=True)
    with BytesIO(archive.raw.data) as b, ZipFile(b) as z:
        z.extractall(DATA)
    # Reassign now
    pdbs = list((DATA / "aligned").glob("**/*_bound.pdb"))
    assert pdbs, "Couldn't find downloaded PDB structures!"
    
print(f"# structures: {len(pdbs)}")


from Bio.PDB import *

# If already split, they should be here:
pdbs_by_chain = list((DATA / "aligned").glob("**/*_bound_chain*[!wH].pdb"))

if not pdbs_by_chain:

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
    print (nr)
    
    # Reassign now
    pdbs_by_chain = list((DATA / "aligned").glob("**/*_bound_chain*.pdb"))

print(f"# structures: {len(pdbs_by_chain)}\n {pdbs_by_chain[0:3]}")


# Protonate

from subprocess import call, STDOUT
def reduce(path):
    path = Path(path)
    output = path.parent / f"{path.stem}_wH{path.suffix}"
    with open(output, "wb") as f:
        call(["reduce", path], stdout=f, stderr=STDOUT)
    return output

from tqdm.auto import tqdm

# If already Hs added, they should be here:
pdbs_by_chain_wH = list((DATA / "aligned").glob("**/*_bound_chain*_wH.pdb"))

if not pdbs_by_chain_wH:
    pdbs_by_chain_wH = {reduce(pdb) for pdb in tqdm(pdbs_by_chain)}

pdbs_by_chain_wH=list(pdbs_by_chain_wH)   
print(f"# structures: {len(pdbs_by_chain_wH)}\n {pdbs_by_chain_wH[0:3]}")

from plipify.core import Structure
from plipify.fingerprints import InteractionFingerprint

structures = []
for path in tqdm(pdbs_by_chain_wH):
    structure = Structure.from_pdbfile(str(path), ligand_name="LIG")
    if len(structure.binding_sites) != 1:
        print(f"{path.relative_to(HERE)} contains {len(structure.binding_sites)} binding sites and we want exactly one.")
        continue
    structures.append(structure)

structures.reverse()


import pandas as pd
lengths = pd.DataFrame([((s.identifier), len(s.sequence()), s.sequence()) for s in structures], columns=["identifier", "length", "sequence"])
# Remove entries where the difference sequence length - median sequence length is greater than one standard deviation
print('Sequence length median and std: ',lengths.length.median(),lengths.length.std() )

lengths["gapcount"]=lengths.sequence.str.count('-')
lengths[lengths.gapcount>0]

# delete those with varying sequence length
lengths = lengths[(lengths.length - lengths.length.median()).abs() < lengths.length.std()]
# delete those with gaps
lengths = lengths[lengths.gapcount==0]
filtered_structures = [s for s in structures if s.identifier in set(lengths.identifier.tolist())]
print('pdbs', len(pdbs_by_chain_wH), "-> structures", len(structures), "-> final", len(filtered_structures), "=", len(pdbs_by_chain_wH) - len(filtered_structures), "structures filtered out")

structures = filtered_structures

from importlib import reload
import plipify.fingerprints, plipify.core, plipify.visualization
reload(plipify.fingerprints)
reload(plipify.core)
reload(plipify.visualization)
from plipify.fingerprints import InteractionFingerprint


# count type of structures
structure_name_type_dict={'Mpro-P':0, 'Mpro-x':0, 'Mpro-z':0, 'other':0}
print(structure_name_type_dict.keys())
for s in structures:
    name_code = s.identifier[0:6]
    if name_code in structure_name_type_dict.keys():
        structure_name_type_dict[name_code]+=1
    else:
        structure_name_type_dict['other']+=1
print(structure_name_type_dict)



# Review
fp = InteractionFingerprint().calculate_fingerprint(
        structures[1:], # see comment above excluded Mpro-z structure
        labeled=True, 
        as_dataframe=True, 
        remove_non_interacting_residues=True,
        remove_empty_interaction_types=True,
        ensure_same_sequence=False
    )

if not fp.values.shape[0]:
    raise ValueError("Fingerprint is empty!")

fp_focused = fp[fp.sum(axis=1) > 5]

# Make publication images with PyMol

selected_structure_name = "Mpro-P0157_0B"
for i,s in enumerate(structures):
    if selected_structure_name in s.identifier:
        print(i, s.identifier)
        selected_structure=s
        selected_structure_id=i
        break
for i, pdb in enumerate(pdbs_by_chain_wH):
    if str(selected_structure_name) in str(pdb):
        print (i, pdb)
        selected_structure_pdb=pdb
        break


#methoxy = Structure.from_pdbfile(path="./data/Mpro/aligned/Mpro-P0157_0A/Mpro-P0157_0A_bound_chainA.pdb")
methoxy = selected_structure

from plipify.visualization import fingerprint_writepdb

output_path = Path("./output/")
pdb_ints = fingerprint_writepdb(fingerprint_df=fp, structure=methoxy, output_path=output_path, ligand=True, summed=True)

from plipify.visualization import VisPymol

#spectrum_cols = ["white_green", "white_pink", "white_cyan", "white_red", "white_blue", "white_yellow", "white_orange", "white_purple"]

spectrum_cols = ["white_cbgreen", "white_cbblue", "white_cborange", "white_cbpurple", "white_cbbrown", "white_pink", "white_cbgrey", "white_yellow"]

# supply view from PyMol "get_view" command

mpro_view = "\
     0.777953982,   -0.552311659,   -0.299568236,\
     0.125412196,   -0.330688834,    0.935369492,\
    -0.615677238,   -0.765243232,   -0.187993601,\
     0.000111578,    0.000063993,  -90.047721863,\
   -22.256666183,    3.898747444,   26.637241364,\
    41.739864349,  138.366195679,  -20.000000000 "

for i, interaction in enumerate(pdb_ints):

    v = VisPymol(pdb=pdb_ints[interaction]._path)
    v.set_style()
    v.create_image(surface=True, ligand_col="green", spectrum_col=spectrum_cols[i], view=mpro_view)
    v.render(name=f"{interaction}_pymol_image", save_path="./output/")


for i, interaction in enumerate(pdb_ints):

    v = VisPymol(pdb=pdb_ints[interaction]._path)
    v.set_style()
    v.create_image(surface=True, ligand_col="green", spectrum_col=spectrum_cols[i], show_ligand=False, view=mpro_view)
    v.render(name=f"{interaction}_pymol_image_no_ligand", save_path="./output/")

# make table

from plipify.visualization import fingerprint_table
html_table = fingerprint_table(fingerprint_df=fp_focused, structure=structures[0])

html_file = open("html_table.html", "w")
html_file.write(html_table)
html_file.close()
