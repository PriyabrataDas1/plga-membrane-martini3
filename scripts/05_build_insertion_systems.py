import os
import math

membrane_gro = "systems/membrane/cancer_membrane.gro"
np_gro = "systems/production_np/production_plga_np_10nm_neg1.gro"

systems = {
    "system_A_partial": {
        "out": "systems/production_combined/system_A_partial/system_A_partial_initial.gro",
        "np_shift_x": 3.0,
        "np_shift_y": 3.0,
        "np_shift_z": 12.0,
        "description": "System A partial insertion, approximately 30-50 percent NP insertion"
    },
    "system_B_deep": {
        "out": "systems/production_combined/system_B_deep/system_B_deep_initial.gro",
        "np_shift_x": 3.0,
        "np_shift_y": 3.0,
        "np_shift_z": 8.0,
        "description": "System B deep insertion, NP center near bilayer center"
    }
}

box_x = 20.0
box_y = 20.0
box_z = 30.0

removal_cutoff = 0.35

def read_gro(filename):
    with open(filename, "r") as f:
        lines = f.readlines()
    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]
    return natoms, atoms

def parse_atom(line):
    resnr = int(line[0:5])
    resname = line[5:10].strip()
    atomname = line[10:15].strip()
    fields = line.split()
    x = float(fields[-3])
    y = float(fields[-2])
    z = float(fields[-1])
    return resnr, resname, atomname, x, y, z

def build_system(name, cfg):
    os.makedirs(os.path.dirname(cfg["out"]), exist_ok=True)

    mem_n, mem_atoms = read_gro(membrane_gro)
    np_n, np_atoms = read_gro(np_gro)

    shifted_np = []
    for line in np_atoms:
        resnr, resname, atomname, x, y, z = parse_atom(line)
        shifted_np.append((
            resnr, resname, atomname,
            x + cfg["np_shift_x"],
            y + cfg["np_shift_y"],
            z + cfg["np_shift_z"]
        ))

    filtered_membrane = []
    removed_beads = 0

    for line in mem_atoms:
        resnr, resname, atomname, x, y, z = parse_atom(line)

        remove = False

        for _, _, _, nx, ny, nz in shifted_np:
            d = math.sqrt((x - nx)**2 + (y - ny)**2 + (z - nz)**2)
            if d < removal_cutoff:
                remove = True
                break

        if remove:
            removed_beads += 1
        else:
            filtered_membrane.append((resnr, resname, atomname, x, y, z))

    combined_atoms = []
    atom_id = 1

    for resnr, resname, atomname, x, y, z in filtered_membrane:
        combined_atoms.append(
            f"{resnr:5d}{resname:<5}{atomname:>5}{atom_id:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
        atom_id += 1

    for resnr, resname, atomname, x, y, z in shifted_np:
        combined_atoms.append(
            f"{resnr:5d}{resname:<5}{atomname:>5}{atom_id:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
        atom_id += 1

    with open(cfg["out"], "w") as f:
        f.write(cfg["description"] + "\n")
        f.write(f"{len(combined_atoms):5d}\n")
        for atom in combined_atoms:
            f.write(atom + "\n")
        f.write(f"{box_x:10.5f}{box_y:10.5f}{box_z:10.5f}\n")

    print("========================================")
    print(name)
    print(f"Output: {cfg['out']}")
    print(f"Original membrane beads: {mem_n}")
    print(f"Removed overlapping beads: {removed_beads}")
    print(f"Remaining membrane beads: {len(filtered_membrane)}")
    print(f"NP beads: {np_n}")
    print(f"Total beads: {len(combined_atoms)}")
    print(f"NP z-shift: {cfg['np_shift_z']}")
    print("========================================")

for name, cfg in systems.items():
    build_system(name, cfg)
