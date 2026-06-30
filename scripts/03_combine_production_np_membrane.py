import os

membrane_gro = "systems/membrane/cancer_membrane.gro"
np_gro = "systems/production_np/production_plga_np_10nm_neg1.gro"

outdir = "systems/production_combined"
out_gro = f"{outdir}/production_plga_membrane_initial.gro"

os.makedirs(outdir, exist_ok=True)

box_x = 20.0
box_y = 20.0
box_z = 30.0

np_shift_x = 3.0
np_shift_y = 3.0
np_shift_z = 17.0

def read_gro(filename):
    with open(filename, "r") as f:
        lines = f.readlines()

    title = lines[0].strip()
    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    box = lines[2 + natoms].strip()

    return title, natoms, atoms, box

mem_title, mem_n, mem_atoms, mem_box = read_gro(membrane_gro)
np_title, np_n, np_atoms, np_box = read_gro(np_gro)

combined_atoms = []

for line in mem_atoms:
    combined_atoms.append(line.rstrip())

atom_id = mem_n + 1

for line in np_atoms:
    resnr = int(line[0:5])
    resname = line[5:10]
    atomname = line[10:15]

    x = float(line[20:28]) + np_shift_x
    y = float(line[28:36]) + np_shift_y
    z = float(line[36:44]) + np_shift_z

    combined_atoms.append(
        f"{resnr:5d}{resname:<5}{atomname:>5}{atom_id:5d}"
        f"{x:8.3f}{y:8.3f}{z:8.3f}"
    )

    atom_id += 1

with open(out_gro, "w") as f:
    f.write("Production PLGA nanoparticle above POPC/POPE/CHOL membrane\n")
    f.write(f"{len(combined_atoms):5d}\n")

    for atom in combined_atoms:
        f.write(atom + "\n")

    f.write(f"{box_x:10.5f}{box_y:10.5f}{box_z:10.5f}\n")

print(f"Combined system written to: {out_gro}")
print(f"Membrane beads: {mem_n}")
print(f"Production NP beads: {np_n}")
print(f"Total beads: {len(combined_atoms)}")
print("Production NP shifted above membrane.")
