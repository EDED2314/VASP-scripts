# FUNCTION EXAMPLES...

# EX 0
# EX 0.1
# generateSlabVac(slab.copy(), "O", 2)
# generateSlabVac(large_slab.copy(), "O", 2)

# Ex 0.2 - adjusting for the top layer of atoms instead of bottom layer (WO2)

# slab = slab.copy()
# slab.set_constraint(None)

# positions = slab.get_positions()
# positions[:, 2] += abs(np.min(positions[:, 2])) + 2.02
# positions[:, 1] += abs(np.min(positions[:, 1]))
# positions[:, 1] = np.abs(positions[:, 1])
# slab.set_positions(positions)

# slab = backup_slab.copy()
# fixed_atoms_indices = get_bottom_n_z_layers(slab, 4)
# constraint = FixAtoms(indices=fixed_atoms_indices)
# slab.set_constraint(constraint)
# print(f"Fixed atoms indices: {fixed_atoms_indices}")
# write("POSCAR_modified", slab)


# EX 1
# newSlab = read("CNST_CONTCAR_H_WO3_TEST")
# fileName = add_h(newSlab, h.copy(), -0.6, "O", 0, dis_x=0, dis_y=-0.60)
# print(fileName)
# generateSimulationFolders(fileName, "2_H")

# fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", 0, pos=(1,2))
# print(fileName)
# generateSimulationFolders(fileName, )

# fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", 0, idxs=triangle_2)
# print(fileName)
# generateSimulationFolders(fileName, )


# for i in range(1, 4):
#     fileName = add_h(
#         large_slab.copy(), h.copy(), height_above_slab, "O", i, dis_x=0, dis_y=0
#     )
#     print(fileName)
#     generateSimulationFolders(
#         fileName, "H_x2y2", templateFolderName="templates_W001_x2y2", trailString="LG"
#     )

# fileName = add_h(slab.copy(), h.copy(), height_above_slab, "W", 0, pos=(0, 0))
# generateSimulationFolders(fileName, trailString="")
# print(fileName)


# for i in range(6):
#     fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", i, layer=-2)
#     generateSimulationFolders(fileName, trailString="L2")
#     print(fileName)

# for i in range(3):
#     fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", i)
#     generateSimulationFolders(fileName, trailString="L1")
#     print(fileName)

# EX 2
# for i in range(3):
#     fileName = add_h2o_vacancy(
#         slab.copy(), h2o.copy(), height_above_slab_for_vacancies, "O", i, "O_down"
#     )
#     print(fileName)
#     generateSimulationFolders(fileName)

# for i in range(1, 4):
#     fileName = add_h2o_vacancy(
#         large_slab.copy(), h2o.copy(), height_above_slab_for_vacancies, "O", i, "O_down"
#     )
#     print(fileName)
#     generateSimulationFolders(
#         fileName, "H2O_x2y2", templateFolderName="templates_W001_x2y2", trailString="LG"
#     )


# EX 3
# UNIT cell - 4 N2 vs 1 N2 vs 2 N2 ... this is for 2 N2 in the unit cell
# large = large_slab.copy()
# fileName = add_n2_vacancy(
#     large, n2.copy(), height_above_slab_for_vacancies, "O", 1, "upright", 0
# )
# fileName = add_n2_vacancy(
#     large, n2.copy(), height_above_slab_for_vacancies, "O", 1, "upright", 0
# )
# print(fileName)
# generateSimulationFolders(
#     fileName, "2N2_x2y2", templateFolderName="templates_W001_x2y2"
# )
# replacePOTCARfromHtoN("2N2_x2y2")

# Ex 3.2 N2 on a WO3 normal unit cell.
# for i in range(3):
#     fileName = add_n2_vacancy(
#         slab.copy(), n2.copy(), height_above_slab_for_vacancies, "O", i, "upright", 0
#     )
#     print(fileName)
#     generateSimulationFolders(fileName)
#     replacePOTCARfromHtoN("N2")


# Ex 4 - adsorbents in vacuums (have to after generation modify the unit cell itself 🤷‍♂️)
# generateAdsorbentInVacuum(emptyCell.copy(), h2o, "H2O")
# generateAdsorbentInVacuum(emptyCell.copy(), h, "H")
# generateAdsorbentInVacuum(emptyCell.copy(), n2, "N2")
# generateAdsorbentInVacuum(emptyCell.copy(), h2, "H2")

# EX 5
# slab = read("CNST_CONTCAR_N2_WO3_V_TEST", format="vasp")
# data, dis = calculateDistancesForEachAtomPair(slab.copy(), "N", "W", 1.37, 0.92)
# slab = read("CNST_CONTCAR_H_WO3_TEST", format="vasp")
# data, dis = calculateDistancesForEachAtomPair(slab.copy(), "O", "H", 0.73, 0.53)
# slab = read("CNST_CONTCAR_N2_WO3_V_TEST", format="vasp")
# data, dis = calculateDistancesForEachAtomPair(slab.copy(), "N", "W")
# slab = read("CNST_CONTCAR_H_WO3_TEST", format="vasp")
# data, dis = calculateDistancesForEachAtomPair(slab.copy(), "O", "H")
# print(data)
# print(dis[0])
# print(dis[1])
# print(dis[2])

# EX 6
# # ex 6.2
# energy, _, _, _ = adsorptionEnergy("OSZICAR_H_WO3", "OSZICAR_WO3", "OSZICAR_H2")
# energy, _, _, _ = adsorptionEnergy("OSZICAR_N2_WO3_V", "OSZICAR_WO3_V", "OSZICAR_N2")
# print(energy)


# EX 7 - visualize

# shutil.rmtree("dataa")
# os.mkdir("dataa")
# fig, ax = plt.subplots()

# for x in [0, 45, 90, 135, 180, 225, 270, 315]:
#     for y in [0, 45, 90, 135, 180, 225, 270, 315]:
#         for z in [0, 45, 90, 135, 180, 225, 270, 315]:
# x = 190
# y = 170
# z = 270
# output_file = f"data/slab_{x}x_{y}y_{z}z.png"
# plot_atoms(slab, ax, rotation=f"{x}x,{y}y,{z}z",)
# ax.set_axis_off()
# plt.savefig(output_file, bbox_inches="tight", pad_inches=0.1, dpi=300)
# # plt.cla()
# plt.close(fig)

# print(f"Slab visualization saved to {output_file}")


# ex 8 - full ex 1

# fileName = add_n2_vacancy(
#     slab.copy(), n2.copy(), height_above_slab, "O", 0, "upright", 0
# )
# print(fileName)
# generateSimulationFolders(
#     fileName, "N2_x2y2_V", templateFolderName="templates_W001_x2y2"
# )
# replacePOTCARfromHtoN("N2_x2y2_V")
# generateSlabVac(large_slab, "O", 0)



#TABLE GENERATION FUNCTIONS

# ---------------------------------------------------------------


def generateHStuff(layer: str):
    post = f"POSTOUTPUT/H_{layer}Layer_OSZICAR"
    post_contcar = f"POSTCONTCAR/H_{layer}Layer_CONTCAR"
    key = NAME_LABEL
    df = pd.DataFrame(
        adsorptionEnergiesOfFolder(
            post,
            "OSZICAR_WO3",
            "OSZICAR_H2",
            multi=0.5,
            name_label=key,
            energy_label=ENERGY_LABEL,
        )
    )
    df = df.sort_values(key)
    df = df.set_index(key)
    # df = df.drop("O0")  # didn't converge yet...
    df = df.drop("P0.0")
    df = df.reset_index()

    df, format_dict = addContcarImagesToDf(df, post_contcar, f"H/{layer}Layer", key)

    refKey = addShortestThreeBondLengthsToDf(df, key, "H", "O", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))
    refKey = addShortestThreeBondLengthsToDf(df, key, "H", "W", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))

    df.to_html(
        "data/H_atom_adsorption_energy.html", escape=False, formatters=format_dict
    )
    print(df)
    return df


def generateH2OStuff(mode: int = 1):
    post = "POSTOUTPUT/H2O_OSZICAR"
    post_contcar = "POSTCONTCAR/H2O_CONTCAR"
    key = NAME_LABEL
    if mode == 1:
        df = pd.DataFrame(
            adsorptionEnergiesOfFolder(
                post,
                "OSZICAR_WO3",
                "OSZICAR_H2",
                name_label=key,
                energy_label=ENERGY_LABEL,
            )
        )
    else:
        df = pd.DataFrame(
            adsorptionEnergiesOfFolder(
                post,
                "OSZICAR_WO3_V_O0",
                "OSZICAR_H2O",
                name_label=key,
                energy_label=ENERGY_LABEL,
            )
        )
    df = df.sort_values(key)
    df = df.set_index(key)
    # df = df.drop("V-O2-OD")  # didn't converge :(
    df = df.reset_index()

    df, format_dict = addContcarImagesToDf(df, post_contcar, "H2O", key, override=False)

    refKey = addShortestThreeBondLengthsToDf(df, key, "H", "O", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))
    refKey = addShortestThreeBondLengthsToDf(df, key, "H", "W", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))
    refKey = addShortestThreeBondLengthsToDf(df, key, "O", "W", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))

    df.to_html(
        f"data/H2O_adsorption_energy_{mode}.html", escape=False, formatters=format_dict
    )
    print(df)
    return df


def generateN2Stuff():
    post = "POSTOUTPUT/N2_OSZICAR"
    post_contcar = "POSTCONTCAR/N2_CONTCAR"
    key = NAME_LABEL
    df = pd.DataFrame(
        adsorptionEnergiesOfFolder(
            post,
            "OSZICAR_WO3_V_O0",
            "OSZICAR_N2",
            name_label=key,
            energy_label=ENERGY_LABEL,
        )
    )
    df = df.sort_values(key)
    df = df.set_index(key)
    df = df.reset_index()

    df, format_dict = addContcarImagesToDf(df, post_contcar, "N2", key)

    refKey = addShortestThreeBondLengthsToDf(df, key, "N", "O", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))
    refKey = addShortestThreeBondLengthsToDf(df, key, "N", "W", post_contcar, "CONTCAR")
    df.insert(2, refKey, df.pop(refKey))

    df.to_html("data/N2_adsorption_energy.html", escape=False, formatters=format_dict)
    print(df)
    return df


# h_wo3_df = generateHStuff("1st")
# h2_wo3_df = generateH2OStuff()
# h2ovwo3_df = generateH2OStuff(mode=2)
# n2_v_wo3_df = generateN2Stuff()

# -------------------------------------------



