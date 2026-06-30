import math
import os

outdir = "systems/nanoparticle"
os.makedirs(outdir, exist_ok=True)

diameter = 10.0
radius = diameter / 2.0
spacing = 0.47

gro_file = f"{outdir}/plga_np_10nm_neg1.gro"
itp_file = f"{outdir}/plga_np_10nm_neg1.itp"

beads = []
idx = 1

n = int(radius / spacing) + 2

for i in range(-n, n + 1):
    for j in range(-n, n + 1):
        for k in range(-n, n + 1):
            x = i * spacing
            y = j * spacing
            z = k * spacing
            r = math.sqrt(x*x + y*y + z*z)
            if r <= radius:
                beads.append((idx, x + radius + 2.0, y + radius + 2.0, z + radius + 2.0, r))
                idx += 1

surface_index = max(beads, key=lambda b: b[4])[0]

with open(gro_file, "w") as f:
    f.write("10 nm negatively charged PLGA nanoparticle MARTINI3\n")
    f.write(f"{len(beads):5d}\n")
    for bead_id, x, y, z, r in beads:
        atomname = "QN" if bead_id == surface_index else "PLG"
        f.write(f"{1:5d}{'PLGA':<5}{atomname:>5}{bead_id:5d}{x:8.3f}{y:8.3f}{z:8.3f}\n")
    f.write(f"{diameter + 4.0:10.5f}{diameter + 4.0:10.5f}{diameter + 4.0:10.5f}\n")

bonds = []
bond_cutoff = spacing * 1.15

for a in range(len(beads)):
    id1, x1, y1, z1, _ = beads[a]
    for b in range(a + 1, len(beads)):
        id2, x2, y2, z2, _ = beads[b]
        d = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        if d <= bond_cutoff:
            bonds.append((id1, id2, d))

with open(itp_file, "w") as f:
    f.write("[ moleculetype ]\n")
    f.write("; name  nrexcl\n")
    f.write("PLGA_NP 1\n\n")

    f.write("[ atoms ]\n")
    f.write("; nr  type  resnr  residue  atom  cgnr  charge  mass\n")
    for bead_id, x, y, z, r in beads:
        if bead_id == surface_index:
            f.write(f"{bead_id:5d} Q5 {1:5d} PLGA QN {bead_id:5d} {-1.0:8.3f} 72.0\n")
        else:
            f.write(f"{bead_id:5d} C2 {1:5d} PLGA PLG {bead_id:5d} {0.0:8.3f} 72.0\n")

    f.write("\n[ bonds ]\n")
    f.write("; i  j  funct  length  force\n")
    for id1, id2, d in bonds:
        f.write(f"{id1:5d} {id2:5d} 1 {d:8.3f} 1250\n")

print(f"Generated {gro_file}")
print(f"Generated {itp_file}")
print(f"Number of beads: {len(beads)}")
print(f"Number of bonds: {len(bonds)}")
print("Net charge: -1")
