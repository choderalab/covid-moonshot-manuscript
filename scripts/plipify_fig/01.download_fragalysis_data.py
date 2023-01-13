import os
from io import BytesIO
from pathlib import Path
from zipfile import ZipFile
import requests
from Bio.PDB import PDBParser, PDBIO

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