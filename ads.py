from typing import List
from ase.io import read, write
from ase.build import add_adsorbate
from ase.build import molecule
from ase import Atom, Atoms
from pymatgen.io.vasp import Poscar, Kpoints
from scipy.spatial.distance import euclidean

import os
from distutils.dir_util import copy_tree
import shutil


class bcolors:
    HEADER = "\033[95m"
    OKBLUE = "\033[94m"
    OKCYAN = "\033[96m"
    OKGREEN = "\033[92m"
    WARNING = "\033[93m"
    FAIL = "\033[91m"
    ENDC = "\033[0m"
    BOLD = "\033[1m"
    UNDERLINE = "\033[4m"


slab = read("CNST_CONTCAR_WO3")
emptyCell = read("CNST_CONTCAR_EMPTY")


# Create the H2O molecule
h2o = molecule("H2O")
n2 = molecule("N2")
n = Atoms("N")
h2 = molecule("H2")
h = Atoms("H")


def cleanUp():
    fileNames = os.listdir()
    for name in fileNames:
        if "POSCAR" in name:
            os.remove(name)
        if "KPOINTS" in name:
            os.remove(name)


def replacePOTCARfromHtoN():
    fileNames = os.listdir()
    N_folders = []
    for name in fileNames:
        if os.path.isdir(name):
            if "N" in name:
                N_folders.append(name)

    for folder in N_folders:
        os.chdir(folder)
        simFolders = os.listdir()
        for sim in simFolders:
            current_potcar_path = os.path.join(sim, "POTCAR")
            shutil.copy("../POTCAR", current_potcar_path)
            print(f"Replaced POTCAR in {sim}")

        os.chdir("..")


def generateSimulationFolders(fileName: str):
    # ex:f"POSCAR_H2O_Vac_{symbol}{index}"
    # ex:f"POSCAR_H2_above_{symbol}{index}"
    # ex:f"POSCAR_N2_Vac_{symbol}{index}_"
    tmp = fileName.split("_")
    idc = tmp[3]  # max 2 char
    symbolRemovedorBelow = idc[0]  # max 1 char cuz its O or W
    moleculeAbove = tmp[1]  # max 3 char
    vac = tmp[2] == "Vac"

    orientation = ""  # max 3 char
    if len(tmp) == 5:
        orientation = tmp[4]

    if not os.path.exists(moleculeAbove):
        os.mkdir(moleculeAbove)

    os.chdir(moleculeAbove)

    if not vac:
        folderName = idc
    else:
        folderName = "V-" + idc

    if orientation != "":
        if orientation == "avg":
            folderName = "Avg-" + idc
        else:
            folderName = folderName + "-" + orientation

    if os.path.exists(folderName):
        shutil.rmtree(folderName)
    os.mkdir(folderName)

    os.chdir("..")

    # template directory containing KPOINTS INCAR and POT and job.slurm
    from_directory = "./templates_W001"
    to_directory = f"./{moleculeAbove}/{folderName}"
    copy_tree(from_directory, to_directory)

    os.rename(fileName, f"./{to_directory}/POSCAR")

    os.chdir(moleculeAbove)
    os.chdir(folderName)

    replacementString = "Et2133JB"  # 2 4 2
    if not vac:
        replacementString = f"{idc}{moleculeAbove}{orientation}"  # 2 3 3
    else:
        if moleculeAbove.lower() == "h2o":
            moleculeAbove = "WT"
        replacementString = f"V{idc}{moleculeAbove}{orientation}"  # 1 2 2 3

    if orientation == "avg":
        replacementString = f"A{idc}{moleculeAbove}"  # 1 4 3

    content = ""
    with open("gpu.slurm", "r") as f:
        content = f.read()

    content = content.replace("JOBNAME", replacementString)

    with open("gpu.slurm", "w") as f:
        f.write(content)

    os.chdir("..")
    os.chdir("..")


def genKpoints(fileName: str):
    try:
        poscar = Poscar.from_file(fileName)
        structure = poscar.structure
        kpoints = Kpoints.automatic_density(structure, kppa=1000)
        newName = fileName.replace("POSCAR_", "")
        kpoints.write_file(f"KPOINTS_{newName}")
    except:
        print("there was an error generating the kpoints file for " + fileName)
        return


def find_average_of_symbol(symbol, idxs):
    points = []

    available_atoms = []
    for atom in slab:
        if atom.symbol == symbol:
            available_atoms.append(atom)

    atom_list = []
    z_max = max([atom.position[2] for atom in available_atoms])
    for atom in available_atoms:
        if abs(atom.position[2] - z_max) < 1e-1:
            atom_list.append(atom.position)

    for idx in idxs:
        points.append(atom_list[idx])

    if not points:
        raise ValueError("The list of points is empty")

    total_x = 0
    total_y = 0
    num_points = len(points)

    for point in points:
        total_x += point[0]
        total_y += point[1]

    avg_x = total_x / num_points
    avg_y = total_y / num_points

    return (avg_x, avg_y)


def remove_atom_at_position(slabb, x, y, atom_type):
    atoms_to_remove = []
    for atom in slabb:
        if (
            atom.symbol == atom_type
            and abs(atom.position[0] - x) < 1e-1
            and abs(atom.position[1] - y) < 1e-1
        ):
            atoms_to_remove.append(atom)

    if len(atoms_to_remove) == 0:
        print("Did not find atom")
        return
    atom = atoms_to_remove[-1]
    slabb.pop(atom.index)


def add_adsorbate_custom(
    slab,
    molecule,
    height,
    symbol,
    index,
    displacement_x=0,
    displacement_y=0,
    vacancy=False,
    idxs=[],
    overridePos=None,
):

    # Determine x,y
    override = not overridePos is None
    if override:
        x = overridePos[0]
        y = overridePos[1]
    else:
        if len(idxs) == 0:
            available_atoms = []
            for atom in slab:
                if atom.symbol == symbol:
                    available_atoms.append(atom)

            if len(available_atoms) == 0:
                print(
                    f'{bcolors.FAIL}Requested symbol: "{symbol}", is not found{bcolors.ENDC}'
                )
                raise ValueError

            atom_list = []
            z_max = max([atom.position[2] for atom in available_atoms])
            for atom in available_atoms:
                if abs(atom.position[2] - z_max) < 1e-1:
                    atom_list.append(atom)

            if len(atom_list) <= index:
                print(
                    f"{bcolors.FAIL}Requested atom index is greater than available atoms{bcolors.ENDC}"
                )
                print(f"{bcolors.FAIL}------{bcolors.ENDC}")
                print(f"Symbol: {bcolors.OKGREEN}{symbol}{bcolors.ENDC}")
                print(
                    f"Length of atom list: {bcolors.OKGREEN}{len(atom_list)}{bcolors.ENDC}"
                )
                print(f"Requested index: {bcolors.OKGREEN}{index}{bcolors.ENDC}")
                for atom in atom_list:
                    print(f"{bcolors.OKBLUE}{atom}{bcolors.ENDC}")
                print(f"{bcolors.FAIL}------{bcolors.ENDC}")
                raise IndexError

            x = atom_list[index].position[0]
            y = atom_list[index].position[1]
            if vacancy:
                remove_atom_at_position(slab, x, y, "O")

        else:
            x, y = find_average_of_symbol(symbol, idxs)

    add_adsorbate(
        slab,
        molecule,
        height,
        (
            x + displacement_x,
            y + displacement_y,
        ),
    )


def add_h(slab, h, height, symbol, index, dis_x=0, dis_y=0, idxs=[], pos=None):
    fileName = f"POSCAR_H_above_{symbol}{index}"
    if 3 >= len(idxs) > 0:
        strIdxs = [str(idx) for idx in idxs]
        symIdx = "".join(strIdxs)
        fileName = f"POSCAR_H_above_{symbol}{symIdx}_avg"
    elif len(idxs) > 3:
        print(
            f"{bcolors.FAIL}Can only do average of three atoms' indices{bcolors.ENDC}"
        )
        return
    add_adsorbate_custom(
        slab, h, height, symbol, index, dis_x, dis_y, idxs=idxs, overridePos=pos
    )
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def add_n(slab, n, height, symbol, index, dis_x=0, dis_y=0, pos=None):
    fileName = f"POSCAR_N_above_{symbol}{index}"
    add_adsorbate_custom(slab, n, height, symbol, index, dis_x, dis_y, overridePos=pos)
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def add_h2(slab, h2, height, symbol, index, dis_x=0, dis_y=0, pos=None):
    fileName = f"POSCAR_H2_above_{symbol}{index}"
    add_adsorbate_custom(slab, h2, height, symbol, index, dis_x, dis_y, overridePos=pos)
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


# Function to add H2O in different orientations
def add_h2o_vacancy(
    slab, h2o, height, symbol, index, orientation="H2_down", rotation=0, pos=None
):
    # 4 -> [4] is the orientation,  2-3 characters
    fileName = f"POSCAR_H2O_Vac_{symbol}{index}_"
    # Center the H2O molecule
    h2o.center()

    if orientation == "H2_down":
        fileName += "H2D"
        pass  # default orientation
    elif orientation == "O_down":
        fileName += "OD"
        h2o.rotate(180, "x")
    elif orientation == "H_down":
        fileName += "HD"
        h2o.rotate(90, "x")
        h2o.rotate(rotation, "z")
        if rotation == 0:
            fileName += "U"
        elif rotation == 90:
            fileName += "L"
        elif rotation == 180:
            fileName += "D"
        elif rotation == 270:
            fileName += "R"
        else:
            fileName += "X"

    elif orientation == "coplanar":
        # add different rotations over here... like coplanar how many deg.
        fileName += "C"
        h2o.rotate(90, "y")
        h2o.rotate(rotation, "z")
        if rotation == 0:
            fileName += "L"
        elif rotation == 90:
            fileName += "D"
        elif rotation == 180:
            fileName += "R"
        elif rotation == 270:
            fileName += "U"
        else:
            fileName += "X"

    add_adsorbate_custom(
        slab, h2o, height, symbol, index, vacancy=True, overridePos=pos
    )
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def add_n2_vacancy(
    slab, n2, height, symbol, index, orientation="upright", rotation=0, pos=None
):
    fileName = f"POSCAR_N2_Vac_{symbol}{index}_"
    # Center the N2 molecule
    n2.center()

    if orientation == "upright":
        fileName += "UPR"
        pass
    elif orientation == "coplanar":
        fileName += "C"
        n2.rotate(90, "y")
        n2.rotate(rotation, "z")
        if rotation == 0:
            fileName += "L"
        elif rotation == 90:
            fileName += "D"
        elif rotation == 180:
            fileName += "R"
        elif rotation == 270:
            fileName += "U"
        else:
            fileName += "X"

    add_adsorbate_custom(slab, n2, height, symbol, index, vacancy=True, overridePos=pos)
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def generateAdsorbentInVacuum(empty, molecule_or_atom, symbol: str):
    fileName = f"POSCAR_{symbol}"
    molecule_or_atom.center()

    empty.pop(0)
    empty += molecule_or_atom

    empty.center(vacuum=20.0)

    write(fileName, empty, format="vasp")

    from_directory = "templates_adsorbate"
    to_directory = f"./adsorbates/{symbol}"

    if os.path.exists(to_directory):
        shutil.rmtree(to_directory)
    os.mkdir(to_directory)

    for template in os.listdir(from_directory):
        if f"POTCAR_{symbol.upper()}" == template:
            shutil.copyfile(
                os.path.join(from_directory, template),
                os.path.join(to_directory, f"POTCAR_{symbol.upper()}"),
            )
        if f"INCAR_{symbol.upper()}" == template:
            shutil.copyfile(
                os.path.join(from_directory, template),
                os.path.join(to_directory, f"INCAR_{symbol.upper()}"),
            )
        if template == "KPOINTS":
            shutil.copyfile(
                os.path.join(from_directory, template),
                os.path.join(to_directory, "KPOINTS"),
            )
        if template == "gpu.slurm":
            shutil.copyfile(
                os.path.join(from_directory, template),
                os.path.join(to_directory, "gpu.slurm"),
            )

    os.rename(fileName, os.path.join(to_directory, fileName))


def calculateDistancesForEachAtomPair(slab, symboll):
    datas = []
    atomss = []
    atomNum = 0

    for i in range(len(slab)):
        for k in range(i + 1, len(slab)):
            # print(slab[i], slab[k])
            data = {}

            point1 = slab[i].position
            point2 = slab[k].position

            data["dis"] = euclidean(point1, point2)
            data["sym1"] = slab[i].symbol
            data["sym2"] = slab[k].symbol
            data["idx1"] = slab[i].index
            data["idx2"] = slab[k].index
            if slab[k].symbol == symboll or slab[i].symbol == symboll:
                datas.append(data)

    # final = []
    # for i in range(len(datas)):
    #     for k in range(i + 1, len(datas)):
    #         if not (
    #             datas[i]["dis"] == datas[k]["dis"]
    #             and datas[i]["dis"] == datas[k]["dis"]
    #         ):
    #             final.append(datas[i])

    return datas


triangle_1 = [0, 1, 4]
triangle_2 = [2, 3, 5]

cleanUp()

# EX 1
# fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", 0, pos=(1,2))
# print(fileName)

# fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", 0, idxs=triangle_2)
# print(fileName)
# generateSimulationFolders(fileName)


# EX 2
# fileName = add_h2o_vacancy(
#     slab.copy(), h2o.copy(), height_above_slab, "O", 0, "coplanar", 270
# )
# print(fileName)
# generateSimulationFolders(fileName)

# EX 3
# fileName = add_n2_vacancy(
#     slab.copy(), n2.copy(), height_above_slab, "O", 0, "coplanar", 0
# )
# print(fileName)
# generateSimulationFolders(fileName)
# replacePOTCARfromHtoN()

# Ex 4
# generateAdsorbentInVacuum(emptyCell.copy(), h2o, "H2O")
# generateAdsorbentInVacuum(emptyCell.copy(), h, "H")
# generateAdsorbentInVacuum(emptyCell.copy(), n2, "N2")

# EX 5
# slab = read("backupPSCR", format="vasp")
# print(slab.get_number_of_atoms())
# data = calculateDistancesForEachAtomPair(slab.copy(), "H")
# print(data)
# for _, pair in enumerate(data):
#     if pair["dis"] < 1:
#         print(pair)


# for i in range(3):
#     fileName = add_h(slab.copy(), h.copy(), height_above_slab, "W", i)
#     print(f"----{i}----")
#     print(fileName)
#     generateSimulationFolders(fileName)

# for i in range(6):
#     fileName = add_h(slab.copy(), h.copy(), height_above_slab, "O", i)
#     print(f"----{i}----")
#     print(fileName)
#     generateSimulationFolders(fileName)


print("----done----")
