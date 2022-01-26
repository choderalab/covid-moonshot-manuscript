# PET-UNK-29afea89-2 Figure

## PET-UNK-29afea89-2 with labelled pockets

*NOTE*: The X-ray code for `PET-UNK-29afea89-2` is `P0157`.

The following command was run to generate the figure:

### Central

```
pymol -cq ../../../scripts/mpro_pocket_flexibility/methoxy_benzo_central.py -- ../../../data/MPro/aligned/ "P0157"
```

The figure was saved as `central_methoxy_benzo_P0157.png`

### RMSDs

```
python ../../../scripts/mpro_pocket_flexibility/calc_pockets_rmsd.py -path ../../../data/MPro/aligned/ -metadata_file ../../../data/MPro/metadata.csv
```
This produces the maximum RMSDs for a particular series in a particular pocket (`mpro_pockets_max_flex_fragments.csv`):

|Pocket | Amino         | Ugi       | Quin | Benzo |
| ------------- | --------- | ---- | ----- | ----- |
|P1| Mpro-x0434_0A | Mpro-P0053_0A | Mpro-x10419_0A | Mpro-x12423_0A |
|P1_prime| Mpro-x12321_0A| Mpro-x2776_0A | Mpro-P0008_0A | Mpro-x10871_0A |
|P2| Mpro-x10598_0A| Mpro-x3359_0A | Mpro-x11366_0A | Mpro-x11757_0A |
|P3_4_5| Mpro-x10019_0A| Mpro-x10082_0A| Mpro-x3303_0A| Mpro-x11424_0A |

with corresponding maximum RMSD values (`mpro_pockets_max_flex_rmsds.csv`):

|Pocket | Amino         | Ugi       | Quin | Benzo |
| ------------- | --------- | ---- | ----- | ----- |
|P1| 1.203598 | 1.056053 | 0.943381 | 0.611776 |
|P1_prime| 0.659057| 0.458427 | 0.501178 | 0.469003 |
|P2| 2.806940| 2.553649 | 2.354497 | 2.630179 |
|P3_4_5| 1.931982| 1.740723| 1.545877| 1.654773 |

Histograms of each ligand series in each pocket can be found in `P1_hist.png`, `P1_prime_hist.png`, `P2.png`, `P3_4_5.png`.

### Pocket flexibility

#### P1 flexibility

The following command was run to generate the PyMol images indicating P1 pocket flexibility:

```python ../../../scripts/mpro_pocket_flexibility/p1_flexibility.py -path ../../../data/MPro/aligned/
```

The following files were produced: `p1_flex_amino_Mpro-x0434_0A.png`, `p1_flex_benzo_Mpro-x12423_0A.png`, `p1_flex_quin_Mpro-x10419_0A.png`, `p1_flex_ugi_Mpro-P0053_0A.png`.

#### P1 prime flexibility

The following command was run to generate the PyMol images indicating P1 prime pocket flexibility:

```python ../../../scripts/mpro_pocket_flexibility/p1_prime_flexibility.py -path ../../../data/MPro/aligned/
```
The following files were produced: `p1_prime_flex_amino_Mpro-x12321_0A.png`, `p1_prime_flex_benzo_Mpro-x10871_0A.png`, `p1_prime_flex_quin_Mpro-P0008_0A.png`, `p1_prime_flex_ugi_Mpro-x2776_0A.png`.

#### P2 flexibility

The following command was run to generate the PyMol images indicating P2 pocket flexibility:

```python ../../../scripts/mpro_pocket_flexibility/p2_flexibility.py -path ../../../data/MPro/aligned/
```
The following files were produced: `p2_flex_amino_Mpro-x10598_0A.png`, `p2_flex_benzo_Mpro-x11757_0A.png`, `p2_flex_quin_Mpro-x11366_0A.png`, `p2_flex_ugi_Mpro-x3359_0A.png`.

#### P3, P4, P5 flexibility

The following command was run to generate the PyMol images indicating P3, P4, P5 pocket flexibility:

```python ../../../scripts/mpro_pocket_flexibility/p3_4_5_flexibility.py -path ../../../data/MPro/aligned/
```
The following files were produced: `p345_flex_amino_Mpro-x10019_0A.png`, `p345_flex_benzo_Mpro-x11424_0A.png`, `p345_flex_quin_Mpro-x3303_0A.png`, `p345_flex_ugi_Mpro-x10082_0A.png`.