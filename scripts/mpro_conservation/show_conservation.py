from sys import argv

import pymol
from pymol import cmd, util

# usage: pymol -cq script_name.py -- pdb_file path/to/fragments/

pdb_arg = argv[1:][0]
print(f"Using PDB file: {pdb_arg} for structures")

pdb_name = f"{pdb_arg.split('.')[0]}"

path = argv[1:][1]

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

cmd.set("surface_smooth_edges", "on")

cmd.reset()
cmd.delete("all")

# Load the central methoxy benzopyran
methoxy_benzo_fragment = 'P0157'
cmd.load(path + f'Mpro-{methoxy_benzo_fragment}_0A/' + f'Mpro-{methoxy_benzo_fragment}_0A.sdf', f'methoxy_benzo-{methoxy_benzo_fragment}-ligand')
cmd.load(f"{pdb_arg}", f'methoxy_benzo-{methoxy_benzo_fragment}-protein')

# Sort out ligand colours
cmd.color("cyan", f"methoxy_benzo-{methoxy_benzo_fragment}-ligand")
util.cbac(f'methoxy_benzo-{methoxy_benzo_fragment}-ligand')

# remove waters
cmd.remove("resn HOH")
cmd.deselect()

# Show molecular representation
cmd.hide("all")
# cmd.dss(f"{pdb_name}")
cmd.bg_color("white")
util.cbaw("*-protein")

# colour based on bfactor
cmd.spectrum("b", "magenta_white", "*-protein", minimum=0, maximum=1)

cmd.show("surface", f"*-protein and not hydrogen")
# cmd.show("ribbon", f"*-protein and not hydrogen")

cmd.disable("*-protein")
cmd.enable(f'methoxy_benzo-{methoxy_benzo_fragment}-protein')

cmd.show("sticks", f"methoxy_benzo-{methoxy_benzo_fragment}-ligand and not hydrogen")

cmd.set("surface_mode", 3)

# Set the viewport and view
cmd.viewport(720, 720)
cmd.set_view("\
     0.885667682,   -0.423480660,   -0.190375701,\
    -0.004551230,   -0.417925209,    0.908467293,\
    -0.464280874,   -0.803735554,   -0.372073978,\
     0.000368759,    0.000287153, -224.885543823,\
   -10.193267822,    6.786962986,   27.136997223,\
    34.983333588,  414.772888184,   20.000000000 ")

pymol.finish_launching()

# Create image
cmd.ray(720, 720)
cmd.png("mpro_conservation.png", dpi=300)
