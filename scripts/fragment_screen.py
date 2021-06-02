import os
import re
from sys import argv

import pymol
from pymol import cmd, util

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

# Load reference aminopyridine structure for alignment
cmd.load(
    path + f"Mpro-x10237_0A/" + f"Mpro-x10237_0A_apo-desolv.pdb",
    f"aminopyridines-x10237-protein",
)

# Load N3 protein and use for alignment
cmd.load("6lu7.pdb")
cmd.align("(6lu7 and chain A)", "aminopyridines-x10237-protein")
cmd.select("N3-protein", "6lu7 and (chain A or chain B)")

# Show all fragments
show_fragments = True
if show_fragments:
    sdf_files = []

    # Get all fragments
    for dirpath, _, filenames in os.walk(f"{path}"):
        for name in filenames:
            if name.endswith(".sdf"):
                # Filter out fragments that aren't in the pocket (manual inspection)
                if name in [
                    "Mpro-x2929_0A.sdf",
                    "Mpro-x1119_0A.sdf",
                    "Mpro-x0887_0A.sdf",
                ]:
                    continue
                sdf_files.append(os.path.join(f"{path}", dirpath, name))

    for sdf_file in sdf_files:
        # Load fragment
        try:
            match = re.search("Mpro-(?P<fragment>x\d+)_", sdf_file)
            fragment = match.group("fragment")
            cmd.load(sdf_file, f"fragment-{fragment}")
        except:
            continue
    cmd.select("fragments", "fragment-*")

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

if show_fragments:
    # Show fragments as lines
    cmd.show("lines", "fragments and not hydrogen")

cmd.set("surface_mode", 3)

# Sort transparency for surfaces and catalytic dyad
cmd.set("transparency", 0.15)
cmd.set("transparency", 1, "resi 145 and resn CYS")
cmd.set("transparency", 1, "resi 41 and resn HIS")

# Set the viewport and view
cmd.viewport(720, 720)
cmd.set_view(
    "\
     0.894281626,   -0.385214210,   -0.227714524,\
     0.184281662,   -0.146693677,    0.971858561,\
    -0.407781065,   -0.911087275,   -0.060196333,\
     0.000453793,    0.000144213,  -90.166122437,\
   -21.955720901,    5.030100346,   28.572034836,\
    76.987503052,  103.323318481,   20.000000000 "
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
cmd.png("fragment_screen_binding_pocket.png")
