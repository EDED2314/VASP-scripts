from typing import List
from ase.io import read, write
from ase.build import add_adsorbate
from ase.build import molecule
from ase import Atom, Atoms
from pymatgen.io.vasp import Poscar, Kpoints
from scipy.spatial.distance import euclidean
from ase.constraints import FixAtoms

import os
from distutils.dir_util import copy_tree
import shutil
import pandas as pd

from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

import numpy as np


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


OUTPUT_DIR = "POSTOUTPUT"


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


def replacePOTCARfromHtoN(N_folder):
    os.chdir(N_folder)
    simFolders = os.listdir()
    for sim in simFolders:
        if os.path.isdir(sim):
            current_potcar_path = os.path.join(sim, "POTCAR")
            shutil.copy("N_POTCAR", current_potcar_path)
            print(f"Replaced POTCAR in {sim}")

    os.chdir("..")


def generateSimulationFolders(
    fileName: str,
    customFolderName="",
    jobFileName="gpu.slurm",
    templateFolderName="templates_W001",
):
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

    if len(customFolderName) == 0:
        mainDirectoryName = moleculeAbove
    else:
        mainDirectoryName = customFolderName

    if not os.path.exists(mainDirectoryName):
        os.mkdir(mainDirectoryName)

    os.chdir(mainDirectoryName)

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
    from_directory = f"./{templateFolderName}"
    to_directory = f"./{mainDirectoryName}/{folderName}"
    shutil.copy(from_directory + "/INCAR", to_directory + "/INCAR")
    shutil.copy(from_directory + "/KPOINTS", to_directory + "/KPOINTS")
    shutil.copy(from_directory + "/POTCAR", to_directory + "/POTCAR")
    shutil.copy(from_directory + "/" + jobFileName, to_directory + "/" + jobFileName)

    os.rename(fileName, f"./{to_directory}/POSCAR")

    os.chdir(mainDirectoryName)
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
    with open(jobFileName, "r") as f:
        content = f.read()

    content = content.replace("JOBNAME", replacementString)

    with open(jobFileName, "w") as f:
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


def getSurfaceAtoms(symbol, index):
    available_atoms = []
    for atom in slab:
        if atom.symbol == symbol:
            available_atoms.append(atom)

    if len(available_atoms) == 0:
        print(f'{bcolors.FAIL}Requested symbol: "{symbol}", is not found{bcolors.ENDC}')
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
        print(f"Length of atom list: {bcolors.OKGREEN}{len(atom_list)}{bcolors.ENDC}")
        print(f"Requested index: {bcolors.OKGREEN}{index}{bcolors.ENDC}")
        for atom in atom_list:
            print(f"{bcolors.OKBLUE}{atom}{bcolors.ENDC}")
        print(f"{bcolors.FAIL}------{bcolors.ENDC}")
        raise IndexError

    return atom_list


def get_bottom_two_layers(slab):
    # Get all z positions of the atoms
    z_positions = [atom.position[2] for atom in slab]
    z_positions.sort()

    # Find the unique z positions and identify the bottom two layers
    unique_z = np.unique(z_positions)
    bottom_two_layers_z = unique_z[:2]

    # Get atoms in the bottom two layers
    bottom_two_layers_atoms = []
    for atom in slab:
        if atom.position[2] in bottom_two_layers_z:
            bottom_two_layers_atoms.append(atom.index)

    return bottom_two_layers_atoms


def generateSlabVac(slab, symbol, index):
    atom_list = getSurfaceAtoms(symbol, index)
    print(atom_list)
    x = atom_list[index].position[0]
    y = atom_list[index].position[1]
    remove_atom_at_position(slab, x, y, "O")
    write("POSCAR", slab, format="vasp")
    genKpoints("POSCAR")


def generateSlab(slab):
    write("POSCAR", slab, format="vasp")
    genKpoints("POSCAR")


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
            atom_list = getSurfaceAtoms(symbol, index)

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
    slab,
    h2o,
    height,
    symbol,
    index,
    orientation="H2_down",
    rotation=0,
    pos=None,
    dis_x=0,
    dis_y=0,
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
        slab,
        h2o,
        height,
        symbol,
        index,
        vacancy=True,
        overridePos=pos,
        displacement_x=dis_x,
        displacement_y=dis_y,
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
    # molecule_or_atom.center()
    molecule_or_atom.center(vacuum=5.0)

    # empty.pop(0)
    # empty += molecule_or_atom

    # empty.center(vacuum=20.0)
    # empty.center()

    write(fileName, molecule_or_atom, format="vasp")
    # write(fileName, empty, format="vasp")

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


# POST SIM ANALYSIS
def adsorptionEnergy(
    OSZICAR_BOTH,
    OSZICAR_SURF,
    OSZICAR_ADS,
    customPathBoth="",
    customPathSurf="",
    customPathAds="",
    adsMulti=1,
):
    """
    First param is for the oszicar of the surface and adsorbate sim, like WO3 vacancy plus H atom
    Second param is the oszicar for the surface
    Third is for the adsorbate
    """
    both = []
    bothDirectory = f"{OUTPUT_DIR}/{OSZICAR_BOTH}"
    if customPathBoth != "":
        bothDirectory = f"{customPathBoth}/{OSZICAR_BOTH}"

    with open(bothDirectory) as f:
        both = f.readlines()

    surf = []
    surfDirectory = f"{OUTPUT_DIR}/{OSZICAR_SURF}"
    if customPathSurf != "":
        surfDirectory = f"{customPathSurf}/{OSZICAR_SURF}"

    with open(surfDirectory) as f:
        surf = f.readlines()

    ads = []
    adsDirectory = f"{OUTPUT_DIR}/{OSZICAR_ADS}"
    if customPathAds != "":
        adsDirectory = f"{customPathAds}/{OSZICAR_ADS}"

    with open(adsDirectory) as f:
        ads = f.readlines()

    lastLineBoth = both[-1]
    lastLineSurf = surf[-1]
    lastLineAds = ads[-1]

    both = lastLineBoth.split()
    surf = lastLineSurf.split()
    ads = lastLineAds.split()

    energyBoth = float(both[2])
    energySurf = float(surf[2])
    energyAds = adsMulti * float(ads[2])  # edit later

    return energyBoth - (energySurf + energyAds)


def analyzeOutputOfFolder(
    POST_DIRECTORY,
    OSZICAR_SURF,
    OSZICAR_ADS,
    multi=1,
    name_label="name",
    energy_label="energy",
):
    datas = []
    for postFile in os.listdir(POST_DIRECTORY):
        data = {}
        data[name_label] = postFile.replace("OSZICAR_", "")
        data[energy_label] = adsorptionEnergy(
            postFile,
            OSZICAR_SURF,
            OSZICAR_ADS,
            customPathBoth=POST_DIRECTORY,
            adsMulti=multi,
        )
        datas.append(data)
    return datas


def calculateDistancesForEachAtomPair(slab, symbol1, symbol2, radius1=0.0, radius2=0.0):
    datas = []
    for i in range(len(slab)):
        for k in range(i + 1, len(slab)):
            data = {}

            point1 = slab[i].position
            point2 = slab[k].position

            data["dis"] = euclidean(point1, point2)
            data["sym1"] = slab[i].symbol
            data["sym2"] = slab[k].symbol
            data["idx1"] = slab[i].index
            data["idx2"] = slab[k].index
            if (slab[k].symbol == symbol1 and slab[i].symbol == symbol2) or (
                slab[k].symbol == symbol2 and slab[i].symbol == symbol1
            ):
                datas.append(data)

    dis = []
    for _, pair in enumerate(datas):
        dis.append(pair["dis"] - (radius1 + radius2))
    dis.sort()

    return datas, dis


slab = read("CNST_CONTCAR_WO3")
# large_slab = read("CNST_CONTCAR_WO3_LG")
emptyCell = read("CNST_CONTCAR_EMPTY")
height_above_slab = 2.2
triangle_1 = [0, 1, 4]
triangle_2 = [2, 3, 5]

cleanUp()

# EX 0
# EX 0.1
# generateSlabVac(slab.copy(), "O", 0)

# Ex 0.2
# slab = read("../backupPSCR", format="vasp")
# fixed_atoms_indices = get_bottom_two_layers(slab)
# constraint = FixAtoms(indices=fixed_atoms_indices)
# slab.set_constraint(constraint)
# write("CONTCAR_fixed", slab)
# print(f"Fixed atoms indices: {fixed_atoms_indices}")

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


# EX 2
# fileName = add_h2o_vacancy(slab.copy(), h2o.copy(), height_above_slab, "O", 0, "O_down")
# print(fileName)
# generateSimulationFolders(fileName)

# EX 3

# # in reality it is no longer O0 but O smth else cuz its a bigger slab
# fileName = add_n2_vacancy(
#     slab.copy(), n2.copy(), height_above_slab, "O", 0, "upright", 0
# )
# print(fileName)
# generateSimulationFolders(
#     fileName, "N2_WO3_x2y2_V", templateFolderName="templates_W001_x2y2"
# )
# replacePOTCARfromHtoN("N2_WO3_x2y2_V")

# Ex 4t
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

# ex 6.1 - post analysis
H_post = "POSTOUTPUT/H_POST"
H_post_contcar = "POSTCONTCAR/H_CONTCAR"
key = "Orientation/Location Molecule Takes"
df = pd.DataFrame(
    analyzeOutputOfFolder(
        H_post,
        "OSZICAR_WO3",
        "OSZICAR_H2",
        multi=0.5,
        name_label=key,
        energy_label="Adsorption Energy (eV)",
    )
)
df = df.sort_values(key)
# df.to_csv("data/adsorption_energy.csv")


def plotThenSaveAtoms(slab, x, y, z, ax, output_file):
    plot_atoms(slab, ax, rotation=f"{x}x,{y}y,{z}z")
    ax.set_axis_off()
    plt.savefig(output_file, bbox_inches="tight", pad_inches=0.1, dpi=300)
    plt.cla()


images = []
if not os.path.exists("images"):
    os.mkdir("images")

fig, ax = plt.subplots()
for fileName in os.listdir(H_post_contcar):
    second_name = fileName.split("_")[1]
    first_name = H_post_contcar.split("/")[1]

    slab = read(f"{H_post_contcar}/{fileName}")

    if not os.path.exists(f"images/{first_name}"):
        os.mkdir(f"images/{first_name}")

    if os.path.exists(f"images/{first_name}/{second_name}"):
        shutil.rmtree(f"images/{first_name}/{second_name}")

    os.mkdir(f"images/{first_name}/{second_name}")

    plotThenSaveAtoms(
        slab,
        135,
        90,
        225,
        ax,
        f"images/{first_name}/{second_name}/slab_135x_90y_225z.png",
    )
    plotThenSaveAtoms(
        slab,
        180,
        180,
        45,
        ax,
        f"images/{first_name}/{second_name}/slab_180x_180y_45z.png",
    )
    plotThenSaveAtoms(
        slab,
        225,
        225,
        35,
        ax,
        f"images/{first_name}/{second_name}/slab_225x_225y_35z.png",
    )

    images.append(
        os.path.abspath(f"images/{first_name}/{second_name}/slab_135x_90y_225z.png")
    )
    images.append(
        os.path.abspath(f"images/{first_name}/{second_name}/slab_180x_180y_45z.png")
    )
    images.append(
        os.path.abspath(f"images/{first_name}/{second_name}/slab_225x_225y_35z.png")
    )

plt.close(fig)


# # your images
# images1 = [
#     "https://vignette.wikia.nocookie.net/2007scape/images/7/7a/Mage%27s_book_detail.png/revision/latest?cb=20180310083825",
#     "https://i.pinimg.com/originals/d9/5c/9b/d95c9ba809aa9dd4cb519a225af40f2b.png",
# ]


# images2 = [
#     "https://static3.srcdn.com/wordpress/wp-content/uploads/2020/07/Quidditch.jpg?q=50&fit=crop&w=960&h=500&dpr=1.5",
#     "https://specials-images.forbesimg.com/imageserve/5e160edc9318b800069388e8/960x0.jpg?fit=scale",
# ]

# df["imageUrls"] = images1
# df["otherImageUrls"] = images2


# # convert your links to html tags
# def path_to_image_html(path):
#     return '<img src="' + path + '" width="60" >'


# pd.set_option("display.max_colwidth", None)

# image_cols = [
#     "imageUrls",
#     "otherImageUrls",
# ]  # <- define which columns will be used to convert to html

# # Create the dictionariy to be passed as formatters
# format_dict = {}
# for image_col in image_cols:
#     format_dict[image_col] = path_to_image_html


df.to_html("data/H_atom_adsorption_energy.html")
print(df)


# # ex 6.2
# energy = adsorptionEnergy("OSZICAR_H_WO3", "OSZICAR_WO3", "OSZICAR_H2")
# energy = adsorptionEnergy("OSZICAR_N2_WO3_V", "OSZICAR_WO3_V", "OSZICAR_N2")
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

# EX8

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

# ex 9 - full ex 1

# fileName = add_n2_vacancy(
#     slab.copy(), n2.copy(), height_above_slab, "O", 0, "upright", 0
# )
# print(fileName)
# generateSimulationFolders(
#     fileName, "N2_WO3_x2y2_V", templateFolderName="templates_W001_x2y2"
# )
# replacePOTCARfromHtoN("N2_WO3_x2y2_V")
# generateSlabVac(large_slab, "O", 0)

print("----done----")
