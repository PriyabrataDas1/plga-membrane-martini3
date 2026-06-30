import math
import os

OUTDIR = "systems/production_np"
os.makedirs(OUTDIR, exist_ok=True)

DIAMETER_NM = 10.0
RADIUS_NM = DIAMETER_NM / 2.0
BEAD_SPACING_NM = 0.47

GRO_FILE = f"{OUTDIR}/production_plga_np_10nm_neg1.gro"
ITP_FILE = f"{OUTDIR}/production_plga_np_10nm_neg1.itp"

beads = []
idx = 1

n = int(RADIUS_NM / BEAD_SPACING_NM) + 2

for i in range(-n, n + 1):
    for j in range(-n, n + 1):
        for k in range(-n, n + 1):
            x = i * BEAD_SPACING_NM
            y = j * BEAD_SPACING_NM
            z = k * BEAD_SPACING_NM
            r = math.sqrt(x*x + y*y + z*z)

            if r <= RADIUS_NM:
                beads.append((idx, x + 7.0, y + 7.0, z + 7.0, r))
                idx += 1

# Assign one surface bead as -1 charge
surface_index = max(beads, key=lambda b: b[4])[0]

with open(GRO_FILE, "w") as f:
    f.write("Production MARTINI 3 PLGA-like NP, 10 nm, net -1\n")
    f.write(f"{len(beads):5d}\n")

    for bead_id, x, y, z, r in beads:
        atomname = "QN" if bead_id == surface_index else "PLG"
        f.write(
            f"{1:5d}{'PLGA':<5}{atomname:>5}{bead_id:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}\n"
        )

    f.write(f"{14.0:10.5f}{14.0:10.5f}{14.0:10.5f}\n")

bonds = []
bond_cutoff = BEAD_SPACING_NM * 1.15

for a in range(len(beads)):
    id1, x1, y1, z1, _ = beads[a]

    for b in range(a + 1, len(beads)):
        id2, x2, y2, z2, _ = beads[b]
        d = math.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

        if d <= bond_cutoff:
            bonds.append((id1, id2, d))

with open(ITP_FILE, "w") as f:
    f.write("; Production PLGA-like MARTINI 3 nanoparticle model\n")
    f.write("; Diameter: 10 nm\n")
    f.write("; Net charge: -1\n")
    f.write("; Neutral core beads: C2\n")
    f.write("; One charged surface bead: Q5, charge -1\n\n")

    f.write("[ moleculetype ]\n")
    f.write("; name    nrexcl\n")
    f.write("PLGA_NP   1\n\n")

    f.write("[ atoms ]\n")
    f.write("; nr  type  resnr  residue  atom  cgnr  charge  mass\n")

    for bead_id, x, y, z, r in beads:
        if bead_id == surface_index:
            f.write(
                f"{bead_id:5d} Q5 {1:5d} PLGA QN "
                f"{bead_id:5d} {-1.0:8.3f} 72.0\n"
            )
        else:
            f.write(
                f"{bead_id:5d} C2 {1:5d} PLGA PLG "
                f"{bead_id:5d} {0.0:8.3f} 72.0\n"
            )

    f.write("\n[ bonds ]\n")
    f.write("; i  j  funct  length  force\n")

    for id1, id2, d in bonds:
        f.write(f"{id1:5d} {id2:5d} 1 {d:8.3f} 1250\n")

print(f"Generated: {GRO_FILE}")
print(f"Generated: {ITP_FILE}")
print(f"Number of beads: {len(beads)}")
print(f"Number of bonds: {len(bonds)}")
print("Net charge: -1")
