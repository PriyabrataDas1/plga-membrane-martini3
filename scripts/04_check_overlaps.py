import math

systems = {
    "surface_control": "systems/production_combined/production_plga_membrane_initial.gro",
    "system_A_partial": "systems/production_combined/system_A_partial/system_A_partial_initial.gro",
    "system_B_deep": "systems/production_combined/system_B_deep/system_B_deep_initial.gro",
}

cutoff_nm = 0.30

def check_overlaps(name, gro_file):
    with open(gro_file, "r") as f:
        lines = f.readlines()

    natoms = int(lines[1].strip())
    atoms = lines[2:2 + natoms]

    membrane_atoms = []
    np_atoms = []

    for line in atoms:
        resname = line[5:10].strip()
        fields = line.split()

        x = float(fields[-3])
        y = float(fields[-2])
        z = float(fields[-1])

        if resname == "PLGA":
            np_atoms.append((x, y, z))
        else:
            membrane_atoms.append((x, y, z))

    min_dist = 999.0
    overlaps = 0

    for nx, ny, nz in np_atoms:
        for mx, my, mz in membrane_atoms:
            d = math.sqrt((nx - mx)**2 + (ny - my)**2 + (nz - mz)**2)

            if d < min_dist:
                min_dist = d

            if d < cutoff_nm:
                overlaps += 1

    print("========================================")
    print(name)
    print(f"File: {gro_file}")
    print(f"NP beads: {len(np_atoms)}")
    print(f"Membrane/solvent beads: {len(membrane_atoms)}")
    print(f"Minimum NP-system distance: {min_dist:.3f} nm")
    print(f"Overlap cutoff: {cutoff_nm:.3f} nm")
    print(f"Number of overlaps: {overlaps}")

    if overlaps == 0:
        print("Status: PASS")
    else:
        print("Status: FAIL")
    print("========================================")

for name, gro_file in systems.items():
    check_overlaps(name, gro_file)
