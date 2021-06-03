from sys import argv

import pymol
from pymol import cmd, util

# usage: pymol -cq script_name.py -- path/to/fragments/

pdb_arg = argv[1:][0]
print(f"Using PDB file: {pdb_arg} for structures")

pdb_name = f"{pdb_arg.split('.')[0]}"

# Set some nice CB friendly colours
cmd.set_color("cb_orange", [0.96, 0.41, 0.23])
cmd.set_color("cb_light_blue", [0.52, 0.75, 0.98])

# Add filter (ambient) and other misc settings
cmd.set("ambient", 0.4)
cmd.set("ambient_occlusion_mode", 1)
cmd.set("ambient_occlusion_scale", 15)
cmd.bg_colour("white")
cmd.set("antialias", 2)
cmd.set("ortho", 1)
cmd.set("ray_trace_mode", 0)

cmd.reset()
cmd.delete("all")

# N3 ligand
cmd.load(f"{pdb_arg}")
cmd.select("N3-ligand", f"{pdb_name} and chain C")
cmd.select("N3-protein", f"{pdb_name} and (chain A or chain B)")
cmd.color("lightpink", "N3-ligand")
util.cnc("N3-ligand")

# remove waters
cmd.remove("resn HOH")
cmd.deselect()

# Show molecular representation
cmd.hide("all")
cmd.dss(f"{pdb_name}")
cmd.bg_color("white")
util.cbaw("*-protein")

# colour based on bfactor
cmd.spectrum("b", "white_green", "*-protein", minimum=0, maximum=1)

# TODO: simplify this selection string
cmd.show(
    "sticks",
    f"*-protein and not hydrogen and not name N and not name C and not name O and not resn 02J and not resn PJE and not resn 010 and not resn PRO and b>0",
)
cmd.show(
    "sticks",
    f"*-protein and not hydrogen and resn PRO and not name C and not name O and b>0",
)
# util.cnc(f"*-protein")

# let beta sheets follow correct path so sticks match up
cmd.set("cartoon_flat_sheets", 0)

cmd.show("cartoon", f"*-protein and not hydrogen")

cmd.show("surface", f"*-protein and not hydrogen")
# cmd.show("ribbon", f"*-protein and not hydrogen")

cmd.disable("*-protein")
cmd.enable("N3-protein")

# cmd.set("surface_color", "white")

molecule = "N3"
cmd.show("sticks", f"{molecule}-ligand and not hydrogen")

cmd.set("surface_mode", 3)

# Sort transparency for surfaces and catalytic dyad
cmd.set("transparency", 0.75)
# cmd.set("transparency", 1, "resi 145 and resn CYS")
# cmd.set("transparency", 1, "resi 41 and resn HIS")

# Set the viewport and view
cmd.viewport(720, 720)
cmd.set_view(
    "\
    -0.201331496,    0.453398257,    0.868266165,\
    -0.455968738,    0.741142511,   -0.492741913,\
    -0.866920233,   -0.495107442,    0.057519738,\
     0.000498764,    0.000351751,  -93.984596252,\
   -10.932416916,   11.968799591,   68.159446716,\
    80.826408386,  107.162223816,   20.000000000 "
)

# Label MPro pockets
cmd.set("label_shadow_mode", 2)

pymol.finish_launching()

# Create image
cmd.ray(720, 720)
cmd.png("mpro_pocket_conservation.png")
