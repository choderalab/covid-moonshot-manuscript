import pymol
from pymol import cmd, util
from sys import argv

# usage: pymol -cq script_name.py -- path/to/fragments/

path = argv[1:][0]
print(f"Using path: {path} for structures")

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
amino_p345 = "Mpro-x10019_0A"
ugi_p345 = "Mpro-x10082_0A"
quin_p345 = "Mpro-x3303_0A"
benzo_p345 = "Mpro-x11424_0A"

fragments = [amino_p345, ugi_p345, quin_p345, benzo_p345]

# Load each fragment: ligand + apo protein
for i, fragment in enumerate(fragments):
    cmd.load(path + f"{fragment}/" + f"{fragment}.sdf", f"p345-{i}-{fragment}-ligand")
    cmd.load(
        path + f"{fragment}/" + f"{fragment}_apo-desolv.pdb",
        f"p345-{i}-{fragment}-protein",
    )

# Sort out colours for each series
cmd.color("wheat", f"p345-0-*")  # amino
cmd.color("palegreen", f"p345-1-*")  # ugi
cmd.color("violet", f"p345-2-*")  # quin
cmd.color("deepolive", f"p345-3-*")  # benzo

# Retain non-carbon default colours
util.cnc(f"p345-*")

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

# P3-5
cmd.select("p3_5", "resi 189+190+191+192+168+167+166+165+192")
cmd.set("surface_color", "palecyan", "p3_5")
cmd.show("surface", "p3_5")
cmd.pseudoatom("p3_5_label", "p3_5")
cmd.set("label_color", "palecyan", "p3_5_label")
cmd.set("label_size", -0.8, "p3_5_label")
cmd.set("label_font_id", 7, "p3_5_label")
### hide psuedoatom
cmd.hide("everything", "p3_5_label")
cmd.show("label", "p3_5_label")

pymol.finish_launching()

# Create images
for i, fragment in zip([0, 1, 2, 3], fragments):
    fragment_key = {0: "amino", 1:"ugi", 2: "quin", 3: "benzo"}
    print(fragment_key[i], fragment)

    cmd.enable(f"p345-{i}-{fragment}-ligand")
    cmd.enable(f"p345-{i}-{fragment}-protein")

    cmd.ray(720, 720)
    cmd.png(f"./p345_flex_{fragment_key[i]}_{fragment}.png", dpi=300)

    cmd.disable(f"p345-{i}-{fragment}-ligand")
    cmd.disable(f"p345-{i}-{fragment}-protein")
