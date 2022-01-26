# MPro conservation Figures

## Creating PDB with conservation 

Before any figures were created an alignment file was used to add conservation information to the `Mpro-P0157_0A_apo-desolv.pdb` structure (obtained from the `MPro` directory in `aligned`. This was achieved via the following command:

```
python ../../scripts/mpro_conservation/add_conservation.py -pdb_file ../../data/mpro_conservation_data/Mpro-P0157_0A_apo-desolv.pdb -alignment_file ../../data/mpro_conservation_data/full_identity_numbers_no_gaps.csv
```

This produced the file `Mpro-P0157_0A_apo-desolv_conservation.pdb`

## Create conservation figures

### Monomer

```
pymol -cq ../../scripts/mpro_conservation/show_conservation.py -- ./Mpro-P0157_0A_apo-desolv_conservation.pdb ../../data/Mpro/aligned/
```

The figure was saved as `mpro_conservation.png`.

### Dimer

```
pymol -cq ../../scripts/mpro_conservation/show_conservation_w_dimer.py -- ./Mpro-P0157_0A_apo-desolv_conservation.pdb ../../data/Mpro/crystallographic/Mpro-P0157.pdb ../../data/MPro/aligned/
```

The figure was saved as `mpro_dimer_conservation.png`.

### Zoomed on active site monomer

```
pymol -cq ../../scripts/mpro_conservation/show_conservation_zoom.py -- ./Mpro-P0157_0A_apo-desolv_conservation.pdb ../../data/MPro/aligned/
```

The figure was saved as `mpro_pocket_conservation.png`.

### Vertical colour bar

Since the color bar cannot be produced in PyMol a custom script can be used to generate the colour bar:

```
python ../../scripts/mpro_conservation/make_cbar.py
```

The final figure was saved as `conervation_cbar.py`.
