import pymol
from pymol import cmd, util
from sys import argv

# usage: pymol -cq script_name.py -- path/to/fragments/

path = argv[1:][0]
print(f"Using path: {path} for structures")

# Set some nice CB friendly colours
cmd.set_color("cb_orange", [0.96, 0.41, 0.23])
cmd.set_color("cb_purple", [0.66, 0.35, 0.63])
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
cmd.load("6lu7.pdb")
cmd.select("N3-ligand", "6lu7 and chain C")
cmd.select("N3-protein", "6lu7 and (chain A or chain B)")
cmd.color("lightpink", "N3-ligand")
util.cnc("N3-ligand")

# remove waters
cmd.remove("resn HOH")
cmd.deselect()

# Show molecular representation
cmd.hide("all")
cmd.dss("6lu7")
cmd.bg_color("white")
util.cbaw("*-protein")

cmd.show("sticks", f"*-protein and not hydrogen")
cmd.show("surface", f"*-protein and not hydrogen")
cmd.disable("*-protein")
cmd.enable("N3-protein")

cmd.set("surface_color", "white")

molecule = "N3"
cmd.show("sticks", f"{molecule}-ligand and not hydrogen")

cmd.set("surface_mode", 3)

# Sort transparency for surfaces and catalytic dyad
cmd.set("transparency", 0.15)
cmd.set("transparency", 1, "resi 145 and resn CYS")
cmd.set("transparency", 1, "resi 41 and resn HIS")

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

cmd.show("sticks", "not hydrogen and *-ligand")

cmd.hide("sticks", "hydrogen")

# Label MPro pockets
cmd.set("label_shadow_mode", 2)

# P1 prime
cmd.select("p1_prime", "resi 25+26+27")
cmd.set("surface_color", "cb_orange", "p1_prime")
cmd.show("surface", "p1_prime")
cmd.pseudoatom("p1_prime_label", "p1_prime")
cmd.set("label_color", "cb_orange", "p1_prime_label")
cmd.set("label_size", -0.8, "p1_prime_label")
cmd.set("label_font_id", 7, "p1_prime_label")
# hide psuedoatom
cmd.hide("everything", "p1_prime_label")
cmd.show("label", "p1_prime_label")

# P1
cmd.select("p1", "resi 142+141+140+172+163+143+144")
cmd.set("surface_color", "cb_yellow", "p1")
cmd.show("surface", "p1")
cmd.pseudoatom("p1_label", "p1")
cmd.set("label_color", "cb_yellow", "p1_label")
cmd.set("label_size", -0.8, "p1_label")
cmd.set("label_font_id", 7, "p1_label")
# hide psuedoatom
cmd.hide("everything", "p1_label")
cmd.show("label", "p1_label")

# P2
cmd.select("p2", "resi 41+49+54")
cmd.set("surface_color", "cb_light_blue", "p2")
cmd.show("surface", "p2")
cmd.pseudoatom("p2_label", "p2")
cmd.set("label_color", "cb_light_blue", "p2_label")
cmd.set("label_size", -0.8, "p2_label")
cmd.set("label_font_id", 7, "p2_label")
# hide psuedoatom
cmd.hide("everything", "p2_label")
cmd.show("label", "p2_label")

# P3-5
cmd.select("p3_5", "resi 189+190+191+192+168+167+166+165+192")
cmd.set("surface_color", "palecyan", "p3_5")
cmd.show("surface", "p3_5")
cmd.pseudoatom("p3_5_label", "p3_5")
cmd.set("label_color", "palecyan", "p3_5_label")
cmd.set("label_size", -0.8, "p3_5_label")
cmd.set("label_font_id", 7, "p3_5_label")
# hide psuedoatom
cmd.hide("everything", "p3_5_label")
cmd.show("label", "p3_5_label")

pymol.finish_launching()

# Create image
cmd.ray(720, 720)
cmd.png("n3_ligand_bound.png", dpi=300)
