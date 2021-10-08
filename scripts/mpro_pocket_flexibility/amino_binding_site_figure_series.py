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

# 3-aminopyridines
# Determined manually from the commented list below (to reduce crowding)
fragments = [
  "Mpro-x11271_0A",
  "Mpro-x11271_0A",
  "Mpro-x11041_0A",
  "Mpro-x10338_0A",
  "Mpro-x10679_0A",
  "Mpro-x10566_0A",
  "Mpro-x11809_0A",
  "Mpro-x0434_0A",
  "Mpro-x10789_0A",
  "Mpro-x11317_0A",
  "Mpro-x10236_0A",
  "Mpro-x10610_0A",
  "Mpro-x0107_0A",
  "Mpro-x10889_1A",
  "Mpro-x10019_0A",
  "Mpro-x10387_0A",
  "Mpro-x2964_0A",
  "Mpro-x10756_0A",
  "Mpro-x10237_0A",
  ]


# fragments = [
#     "Mpro-x11271_0A",
#     "Mpro-x11041_0A",
#     "Mpro-x10338_0A",
#     "Mpro-x10555_0A",
#     "Mpro-x10679_0A",
#     "Mpro-x10566_0A",
#     "Mpro-x10889_0A",
#     "Mpro-x11809_0A",
#     "Mpro-x0434_0A",
#     "Mpro-x10789_0A",
#     "Mpro-x11317_0A",
#     "Mpro-x10236_0A",
#     "Mpro-x10610_0A",
#     "Mpro-x11723_0A",
#     "Mpro-x11044_0A",
#     "Mpro-x10942_0A",
#     "Mpro-x11473_0A",
#     "Mpro-x2912_0A",
#     "Mpro-x0107_0A",
#     "Mpro-x11233_0A",
#     "Mpro-x10334_0A",
#     "Mpro-x10889_1A",
#     "Mpro-x11186_0A",
#     "Mpro-x10638_0A",
#     "Mpro-x10377_0A",
#     "Mpro-x11485_0A",
#     "Mpro-x11801_0A",
#     "Mpro-x11346_0A",
#     "Mpro-x11562_0A",
#     "Mpro-x0678_0A",
#     "Mpro-x10019_0A",
#     "Mpro-x10324_0A",
#     "Mpro-x2569_0A",
#     "Mpro-x10395_0A",
#     "Mpro-x10473_0A",
#     "Mpro-x11011_0A",
#     "Mpro-x11318_0A",
#     "Mpro-x3366_0A",
#     "Mpro-x10898_0A",
#     "Mpro-x11013_0A",
#     "Mpro-x11225_0A",
#     "Mpro-x10906_0A",
#     "Mpro-x10565_0A",
#     "Mpro-x10476_0A",
#     "Mpro-x2581_0A",
#     "Mpro-x10834_0A",
#     "Mpro-x10387_0A",
#     "Mpro-x10535_0A",
#     "Mpro-x11475_0A",
#     "Mpro-x10645_0A",
#     "Mpro-x10494_0A",
#     "Mpro-x10604_0A",
#     "Mpro-x10787_0A",
#     "Mpro-x2643_0A",
#     "Mpro-x10976_0A",
#     "Mpro-x3298_0A",
#     "Mpro-x10900_0A",
#     "Mpro-x10488_0A",
#     "Mpro-x2646_0A",
#     "Mpro-x10314_0A",
#     "Mpro-x2908_0A",
#     "Mpro-x10888_0A",
#     "Mpro-x2971_0A",
#     "Mpro-x11159_0A",
#     "Mpro-x10723_0A",
#     "Mpro-x10417_0A",
#     "Mpro-x10801_0A",
#     "Mpro-x11231_0A",
#     "Mpro-x10359_0A",
#     "Mpro-x10559_0A",
#     "Mpro-x10423_0A",
#     "Mpro-x11428_0A",
#     "Mpro-x2600_0A",
#     "Mpro-x2964_0A",
#     "Mpro-x10392_0A",
#     "Mpro-x10247_0A",
#     "Mpro-x10201_0A",
#     "Mpro-x10756_0A",
#     "Mpro-x3108_0A",
#     "Mpro-x11368_0A",
#     "Mpro-x10022_0A",
#     "Mpro-x10329_0A",
#     "Mpro-x3080_0A",
#     "Mpro-x10484_0A",
#     "Mpro-x10606_0A",
#     "Mpro-x2572_0A",
#     "Mpro-x11417_0A",
#     "Mpro-x11164_0A",
#     "Mpro-x10248_0A",
#     "Mpro-x2562_0A",
#     "Mpro-x2649_0A",
#     "Mpro-x10996_0A",
#     "Mpro-x10422_0A",
#     "Mpro-x10396_0A",
#     "Mpro-x10237_0A",
#     "Mpro-x11354_0A",
#     "Mpro-x10322_0A",
#     "Mpro-x11372_0A",
#     "Mpro-x11339_0A",
#     "Mpro-x10995_0A",
#     "Mpro-x10733_0A",
#     "Mpro-x11564_0A",
#     "Mpro-x10474_0A",
#     "Mpro-x11426_0A",
#     "Mpro-x10575_0A",
#     "Mpro-x11488_0A",
#     "Mpro-x10478_0A",
#     "Mpro-x11427_0A",
#     "Mpro-x11743_0A",
#     "Mpro-x10598_0A",
#     "Mpro-x2608_0A",
#     "Mpro-x10856_0A",
#     "Mpro-x10421_0A",
#     "Mpro-x11025_0A",
#     "Mpro-x10626_0A",
#     "Mpro-x10800_0A",
#     "Mpro-x10371_0A",
#     "Mpro-x10178_0A",
#     "Mpro-x10327_0A",
#     "Mpro-x11764_0A",
#     "Mpro-x11532_0A",
#     "Mpro-x12321_0A",
#     "Mpro-x12300_0A",
#     "Mpro-x11507_0A",
#     "Mpro-x11557_0A",
#     "Mpro-x11708_0A",
#     "Mpro-x11641_0A",
#     "Mpro-x11540_0A",
#     "Mpro-x11560_0A",
#     "Mpro-x11513_0A",
#     "Mpro-x11543_0A",
#     "Mpro-x12682_0A",
#     "Mpro-x12740_0A",
#     "Mpro-x12719_0A",
#     "Mpro-P0045_0A",
# ]

# Load each fragment: ligand + apo protein
for fragment in fragments:
    cmd.load(
        path + f"{fragment}/" + f"{fragment}.sdf", f"aminopyridines-{fragment}-ligand"
    )
    cmd.load(
        path + f"{fragment}/" + f"{fragment}_apo-desolv.pdb",
        f"aminopyridines-{fragment}-protein",
    )
cmd.color("wheat", f"aminopyridines-*")
util.cnc(f"aminopyridines-*")

# remove waters
cmd.remove("resn HOH")
cmd.deselect()

# Show molecular representation
cmd.hide("all")
cmd.bg_color("white")
util.cbaw("*-protein")

cmd.show("sticks", f"*-protein and not hydrogen")
cmd.show("surface", f"aminopyridines-Mpro-x11271_0A-protein and not hydrogen")
cmd.disable("*-protein")
cmd.enable("aminopyridines-Mpro-x11271_0A-protein and not hydrogen")

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
cmd.png("./amino_lead_series.png", dpi=300)
