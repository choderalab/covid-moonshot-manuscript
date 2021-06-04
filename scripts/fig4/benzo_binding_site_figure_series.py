import pymol
from pymol import cmd, util
from sys import argv

# usage: pymol -cq script_name.py -- path/to/fragments/

path = argv[1:][0]
print(f"Using path: {path} for structures")

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

fragments = [
    "Mpro-x10820_0A",
    "Mpro-x10876_0A",
    "Mpro-x10870_0A",
    "Mpro-x11454_0A",
    "Mpro-x11493_0A",
    "Mpro-x11424_0A",
    "Mpro-x11432_0A",
    "Mpro-x10871_0A",
    "Mpro-x10466_0A",
    "Mpro-x11798_0A",
    "Mpro-x11797_0A",
    "Mpro-x12423_0A",
    "Mpro-x12177_0A",
    "Mpro-x11789_0A",
    "Mpro-x12204_0A",
    "Mpro-x12026_0A",
    "Mpro-x12419_0A",
    "Mpro-x12143_0A",
    "Mpro-x11757_0A",
    "Mpro-x12136_0A",
]

# Load each fragment: ligand + apo protein
for fragment in fragments:
    cmd.load(
        path + f"{fragment}/" + f"{fragment}.sdf", f"benzotriazoles-{fragment}-ligand"
    )
    cmd.load(
        path + f"{fragment}/" + f"{fragment}_apo-desolv.pdb",
        f"benzotriazoles-{fragment}-protein",
    )
cmd.color("deepolive", f"benzotriazoles-*")
util.cnc(f"benzotriazoles-*")

# remove waters
cmd.remove("resn HOH")
cmd.deselect()

# Show molecular representation
cmd.hide("all")
cmd.bg_color("white")
util.cbaw("*-protein")

cmd.show("sticks", f"*-protein and not hydrogen")
cmd.show("surface", f"benzotriazoles-Mpro-x10820_0A-protein and not hydrogen")
cmd.disable("*-protein")
cmd.enable("benzotriazoles-Mpro-x10820_0A-protein and not hydrogen")

cmd.set("surface_color", "white")

cmd.set("surface_mode", 3)
cmd.set("transparency", 0.15)
cmd.set("transparency", 1, "resi 145 and resn CYS")
cmd.set("transparency", 1, "resi 41 and resn HIS")

cmd.viewport(720, 720)

cmd.set_view(
    "\
     0.889428616,   -0.456408978,    0.024401871,\
    -0.099697880,   -0.141634151,    0.984882414,\
    -0.446053863,   -0.878418088,   -0.171478689,\
     0.000357639,    0.000210542,  -94.736244202,\
   -21.802495956,    4.796146870,   28.488866806,\
    81.561714172,  107.897544861,   20.000000000 "
)

cmd.show("sticks", "not hydrogen and *-ligand")

# Hide all fragments
cmd.deselect()

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
# hide pseudoatom
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
# hide pseudoatom
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
# hide pseudoatom
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
# hide pseudoatom
cmd.hide("everything", "p3_5_label")
cmd.show("label", "p3_5_label")

pymol.finish_launching()

cmd.ray(720, 720)
cmd.png("./benzo_lead_series.png", dpi=300)
