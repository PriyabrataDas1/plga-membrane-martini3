systems = {
    "Surface": "runs/surface_control/npt/surface_control_npt.gro",
    "System_A": "runs/system_A/npt/system_A_npt.gro",
    "System_B": "runs/system_B/npt/system_B_npt.gro",
    "Partial_3500ps": "systems/pulled_systems/system_partial_3500ps.gro",
}

for label, filename in systems.items():

    with open(filename) as f:
        lines = f.readlines()

    natoms = int(lines[1])

    plga_z = []
    lipid_z = []

    for line in lines[2:2 + natoms]:

        resname = line[5:10].strip()
        z = float(line[36:44])

        if resname == "PLGA":
            plga_z.append(z)

        elif resname in ["POPC", "POPE", "CHOL"]:
            lipid_z.append(z)

    np_com = sum(plga_z) / len(plga_z)
    mem_center = sum(lipid_z) / len(lipid_z)
    distance = np_com - mem_center

    print("=" * 50)
    print(label)
    print(f"NP COM z        : {np_com:.3f} nm")
    print(f"Membrane center : {mem_center:.3f} nm")
    print(f"Distance        : {distance:.3f} nm")
