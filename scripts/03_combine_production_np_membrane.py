import os
import math

membrane_gro = "systems/membrane/cancer_membrane.gro"
np_gro = "systems/production_np/production_plga_np_10nm_neg1.gro"

outdir = "systems/production_combined"
out_gro = f"{outdir}/production_plga_membrane_initial.gro"

os.makedirs(outdir, exist_ok=True)

BOX = (20.0, 20.0, 30.0)

np_shift_x = 3.0
np_shift_y = 3.0
np_shift_z = 17.0

water_removal_cutoff = 0.35


def read_gro(filename):
    with open(filename) as f:
        lines = f.readlines()

    natoms = int(lines[1].strip())
    atoms = []

    for line in lines[2:2 + natoms]:
        fields = line.split()
        atoms.append({
            "resnr": int(line[0:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "x": float(fields[-3]),
            "y": float(fields[-2]),
            "z": float(fields[-1]),
        })

    return atoms


def distance(a, b):
    return math.sqrt(
        (a["x"] - b["x"]) ** 2 +
        (a["y"] - b["y"]) ** 2 +
        (a["z"] - b["z"]) ** 2
    )


def write_atom(f, atom, atomnr):
    f.write(
        f"{atom['resnr'] % 100000:5d}"
        f"{atom['resname'][:5]:<5}"
        f"{atom['atomname'][:5]:>5}"
        f"{atomnr % 100000:5d}"
        f"{atom['x']:8.3f}"
        f"{atom['y']:8.3f}"
        f"{atom['z']:8.3f}\n"
    )


membrane = read_gro(membrane_gro)
np_atoms = read_gro(np_gro)

shifted_np = []

for atom in np_atoms:
    shifted_np.append({
        "resnr": atom["resnr"],
        "resname": atom["resname"],
        "atomname": atom["atomname"],
        "x": atom["x"] + np_shift_x,
        "y": atom["y"] + np_shift_y,
        "z": atom["z"] + np_shift_z,
    })

filtered_membrane = []
removed_water = 0

for atom in membrane:
    is_water = atom["resname"] == "W" or atom["atomname"] == "W"

    remove = False

    if is_water:
        for npa in shifted_np:
            if distance(atom, npa) < water_removal_cutoff:
                remove = True
                break

    if remove:
        removed_water += 1
    else:
        filtered_membrane.append(atom)

all_atoms = filtered_membrane + shifted_np

with open(out_gro, "w") as f:
    f.write("Surface-control PLGA NP above POPC/POPE/CHOL membrane; waters removed near NP\n")
    f.write(f"{len(all_atoms):5d}\n")

    atomnr = 1

    for atom in all_atoms:
        write_atom(f, atom, atomnr)
        atomnr += 1

    f.write(f"{BOX[0]:10.5f}{BOX[1]:10.5f}{BOX[2]:10.5f}\n")

print(f"Output: {out_gro}")
print(f"Original membrane beads: {len(membrane)}")
print(f"Waters removed: {removed_water}")
print(f"Remaining membrane/solvent beads: {len(filtered_membrane)}")
print(f"NP beads: {len(shifted_np)}")
print(f"Total atoms: {len(all_atoms)}")
print(f"Expected topology atoms: {len(all_atoms)}")
