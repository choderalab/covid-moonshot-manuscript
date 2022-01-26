# covid-moonshot-manuscript
Scripts, figures, and data relevant to the COVID moonshot manuscript.

# Organisation

`scripts`: contains PyMol and Python scripts used to create protein - ligand images.
In order to run these scripts the aligned structures and `metadata.csv` are required. These can be downloaded via: 

```
wget https://fragalysis.diamond.ac.uk/media/targets/Mpro.zip
```
It is recommened that the `.zip` file is extracted in the `./data` directory of this repository since that is where all scripts in the `./figures` `README.md` files are referenced to.

## Mpro data structure

* Run `unzip MPro.zip` to extract the relevant files.
* The `aligned/` directory contains all aligned MPro structures.
* `metadata.csv` contains relevant information for fragment names etc.

# Examples

* Check the `README.md` appropriate files in the `./figures` directory for example terminal commands used to create each figure.
* A working installation of [PyMol](https://pymol.org/2/) is required.
* It is also recommended to have the following line in your `./pymolrc` to access colour blind friendly colours:
```
run https://github.com/Pymol-Scripts/Pymol-script-repo/raw/master/colorblindfriendly.py
```



