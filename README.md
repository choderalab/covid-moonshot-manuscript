# covid-moonshot-manuscript
Scripts, figures, and data relevant to the COVID moonshot manuscript

# Organisation

`scripts`: contains PyMol and Python scripts used to create protein - ligand images.
 In order to run these scripts the aligned structures and `metadata.csv` are required. These can be downloaded via: 
 ```
 wget https://fragalysis.diamond.ac.uk/media/targets/Mpro.zip
 ```
 
 * Run `unzip MPro.zip` to extract the relevant files.
 * The `aligned/` directory contains all aligned MPro structures.
 * `metadata.csv` contains relevant information for fragment names etc.

# Examples

* To generate an image of the fragment screen at the MPro binding site (Fig 1d) run:
```
pymol -cq scripts/fig1/script_name.py -- path/to/mpro/data/aligned/
```
* To calculate the RMSD of each pocket (Fig 4b) run: 
```
python calc_pockets_rmsd.py -path path/to/fragments/aligned -metadata_file path/to/mpro/data/metadata.csv
```


