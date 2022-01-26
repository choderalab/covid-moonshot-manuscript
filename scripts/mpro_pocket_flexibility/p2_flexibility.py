import pymol
from pymol import cmd, util
from sys import argv

# usage: pymol -cq script_name.py -- path/to/fragments/

path = argv[1:][0]
print(f"Using path: {path} for structures")

# Set some nice CB friendly colours
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

# Specify fragment codes from each series
amino_p2 = "Mpro-x10598_0A"
ugi_p2 = "Mpro-x3359_0A"
quin_p2 = "Mpro-x11366_0A"
benzo_p2 = "Mpro-x11757_0A"

fragments = [amino_p2, ugi_p2, quin_p2, benzo_p2]

# Load each fragment: ligand + apo protein
for i, fragment in enumerate(fragments):
    cmd.load(path + f"{fragment}/" + f"{fragment}.sdf", f"p2-{i}-{fragment}-ligand")
    cmd.load(
        path + f"{fragment}/" + f"{fragment}_apo-desolv.pdb",
        f"p2-{i}-{fragment}-protein",
    )

# Sort out colours for each series
cmd.color("wheat", f"p2-0-*")  # amino
cmd.color("palegreen", f"p2-1-*")  # ugi
cmd.color("violet", f"p2-2-*")  # quin
cmd.color("deepolive", f"p2-3-*")  # benzo

# Retain non-carbon default colours
util.cnc(f"p2-*")

# remove waters
cmd.remove("resn HOH")
cmd.deselect()

# Show molecular representation
cmd.hide("all")
cmd.bg_color("white")
util.cbaw("*-protein")

cmd.show("sticks", f"*-protein and not hydrogen")
cmd.show("surface", f"*-protein and not hydrogen")
cmd.show("sticks", "not hydrogen and *-ligand")

cmd.hide("sticks", "hydrogen")

cmd.disable("*-protein")
cmd.disable("*-ligand")

cmd.set("surface_color", "white")

# Sort transparency for surfaces and catalytic dyad
cmd.set("surface_mode", 3)
cmd.set("transparency", 0.15)
cmd.set("transparency", 1, "resi 145 and resn CYS")
cmd.set("transparency", 1, "resi 41 and resn HIS")

# Set the viewport and view
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

# Label MPro pockets
cmd.set("label_shadow_mode", 2)

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

pymol.finish_launching()

# Create images
print("Creating images for P2...")
for i, fragment in zip([0, 1, 2, 3], fragments):
    fragment_key = {0: "amino", 1:"ugi", 2: "quin", 3: "benzo"}
    print(fragment_key[i], fragment)

    cmd.enable(f"p2-{i}-{fragment}-ligand")
    cmd.enable(f"p2-{i}-{fragment}-protein")

    cmd.ray(720, 720)
    cmd.png(f"./p2_flex_{fragment_key[i]}_{fragment}.png", dpi=300)

    cmd.disable(f"p2-{i}-{fragment}-ligand")
    cmd.disable(f"p2-{i}-{fragment}-protein")
