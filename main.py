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

pd.set_option("display.max_colwidth", None)

from ase.visualize.plot import plot_atoms
import matplotlib.pyplot as plt

import numpy as np

import matplotlib.pyplot as plt
from pymatgen.io.vasp import Vasprun

from energies import *
from constants import *


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
        elif "KPOINTS" in name:
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
    trailString="",
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

    replacementString += trailString

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


def find_average_of_symbol(symbol, idxs, slab, layer):
    points = []

    atom_list = getSurfaceAtoms(symbol, 0, slab, layer=layer)
    atom_list = [atom.position for atom in atom_list]

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


def remove_atom_at_position_on_surface(slabb, x, y, atom_type):
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

    atom_selected = atoms_to_remove[0]
    for atom_want in atoms_to_remove:
        if atom_want.position[2] > atom_selected.position[2]:
            atom_selected = atom_want

    slabb.pop(atom_selected.index)


def getSurfaceAtoms(symbol: str, index: int, slab, layer: int = -1):
    available_atoms = []
    for atom in slab:
        if atom.symbol == symbol:
            available_atoms.append(atom)

    if len(available_atoms) == 0:
        print(f'{bcolors.FAIL}Requested symbol: "{symbol}", is not found{bcolors.ENDC}')
        raise ValueError

    atom_list = []
    z_pos = [atom.position[2] for atom in available_atoms]
    # Round the numbers to the desired precision and then use a set to remove duplicates
    rounded_list = [round(num, 3) for num in z_pos]
    z_pos = list(set(rounded_list))
    z_pos.sort()
    z_wanted = z_pos[layer]
    for atom in available_atoms:
        if abs(atom.position[2] - z_wanted) < 1e-1:
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


def get_bottom_n_z_layers(slab, n: int):
    # Get all z positions of the atoms
    z_positions = [atom.position[2] for atom in slab]
    z_positions.sort()

    # Find the unique z positions and identify the bottom two layers

    unique_z = np.unique(z_positions)
    bottom_two_layers_z = unique_z[:n]

    # Get atoms in the bottom two layers
    atoms = []
    for atom in slab:
        if atom.position[2] in bottom_two_layers_z:
            atoms.append(atom.index)

    return atoms


def generateSlabVac(slab, symbol, index):
    atom_list = getSurfaceAtoms(symbol, index, slab)
    x = atom_list[index].position[0]
    y = atom_list[index].position[1]
    remove_atom_at_position_on_surface(slab, x, y, symbol)
    write("POSCAR", slab, format="vasp")
    genKpoints("POSCAR")


def generateSlab(slab):
    write("POSCAR", slab, format="vasp")
    genKpoints("POSCAR")


def addAdsorbateCustom(
    slab,
    molecule,
    height: float,
    symbol: str,
    index: int,
    displacement_x: float = 0,
    displacement_y: float = 0,
    vacancy=False,
    idxs=[],
    overridePos=None,
    layer: int = -1,
):

    # Determine x,y
    override = not overridePos is None
    if override:
        x = overridePos[0]
        y = overridePos[1]
    else:
        if len(idxs) == 0:
            atom_list = getSurfaceAtoms(symbol, index, slab, layer)

            x = atom_list[index].position[0]
            y = atom_list[index].position[1]
            if vacancy:
                remove_atom_at_position_on_surface(slab, x, y, "O")

        else:
            x, y = find_average_of_symbol(symbol, idxs, slab, layer)

    add_adsorbate(
        slab,
        molecule,
        height,
        (
            x + displacement_x,
            y + displacement_y,
        ),
    )


def add_h(
    slab, h, height, symbol, index, dis_x=0, dis_y=0, idxs=[], pos=None, layer=-1
):
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
    addAdsorbateCustom(
        slab,
        h,
        height,
        symbol,
        index,
        dis_x,
        dis_y,
        idxs=idxs,
        overridePos=pos,
        layer=layer,
    )
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def add_n(slab, n, height, symbol, index, dis_x=0, dis_y=0, pos=None):
    fileName = f"POSCAR_N_above_{symbol}{index}"
    addAdsorbateCustom(slab, n, height, symbol, index, dis_x, dis_y, overridePos=pos)
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def add_h2(slab, h2, height, symbol, index, dis_x=0, dis_y=0, pos=None):
    fileName = f"POSCAR_H2_above_{symbol}{index}"
    addAdsorbateCustom(slab, h2, height, symbol, index, dis_x, dis_y, overridePos=pos)
    write(fileName, slab, format="vasp")
    genKpoints(fileName)
    return fileName


def add_h2o_to_existing_configurations_from_directory():

    return


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

    addAdsorbateCustom(
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
        n2.rotate(180, "x")
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

    addAdsorbateCustom(slab, n2, height, symbol, index, vacancy=True, overridePos=pos)
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
def plotCompleteDOS(listOfVaspRunFilePaths, colorsList, sigma=0.1):
    assert len(listOfVaspRunFilePaths) == len(colorsList)

    for i in range(len(listOfVaspRunFilePaths)):
        path = listOfVaspRunFilePaths[i]
        color = colorsList[i]
        vr = Vasprun(path)

        dos = vr.complete_dos
        dos.densities = dos.get_smeared_densities(sigma)

        x = dos.energies - dos.efermi
        y = dos.get_densities()

        plt.plot(x, y, color=color)

    plt.xlabel("E - Ef (eV)")
    plt.ylabel("DOS (a.u.)")
    plt.show()


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
    bothDirectory = f"{OUTPUT_DIR}/{OSZICAR_BOTH}"
    if customPathBoth != "":
        bothDirectory = f"{customPathBoth}/{OSZICAR_BOTH}"

    surfDirectory = f"{OUTPUT_DIR}/{OSZICAR_SURF}"
    if customPathSurf != "":
        surfDirectory = f"{customPathSurf}/{OSZICAR_SURF}"

    adsDirectory = f"{OUTPUT_DIR}/{OSZICAR_ADS}"
    if customPathAds != "":
        adsDirectory = f"{customPathAds}/{OSZICAR_ADS}"

    energyBoth = readOszicarFileAndGetLastLineEnergy(bothDirectory)
    energySurf = readOszicarFileAndGetLastLineEnergy(surfDirectory)
    energyAds = adsMulti * readOszicarFileAndGetLastLineEnergy(adsDirectory)

    return energyBoth - (energySurf + energyAds), energyBoth, energySurf, energyAds


def adsorptionEnergiesOfFolder(
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
        data[energy_label], _, _, _ = adsorptionEnergy(
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
            if (slab[k].symbol == symbol1 and slab[i].symbol == symbol2) or (
                slab[k].symbol == symbol2 and slab[i].symbol == symbol1
            ):
                data = {}

                point1 = slab[i].position
                point2 = slab[k].position

                data["dis"] = euclidean(point1, point2)
                data["sym1"] = slab[i].symbol
                data["sym2"] = slab[k].symbol
                data["idx1"] = slab[i].index
                data["idx2"] = slab[k].index

                datas.append(data)

    dis = []
    for _, pair in enumerate(datas):
        dis.append(pair["dis"])
    dis.sort()

    return datas, dis


def addContcarImagesToDf(
    df, CONTCAR_DIRECTORY: str, POSCAR_DIRECTORY: str, key: str, override=False
):

    def plotThenSaveAtoms(slab, x, y, z, ax, output_file):
        plot_atoms(slab, ax, rotation=f"{x}x,{y}y,{z}z")
        ax.set_axis_off()
        plt.savefig(output_file, bbox_inches="tight", pad_inches=0.1, dpi=300)
        plt.cla()

    def path_to_image_html(path):
        return '<img src="' + path + '" width="200" >'

    images1 = []
    images2 = []
    images3 = []
    initImages1 = []
    initImages2 = []
    initImages3 = []
    if not os.path.exists("images"):
        os.mkdir("images")

    names = df[key]
    fig, ax = plt.subplots()
    for name in names:
        initPoscar = f"{POSCAR_DIRECTORY}/{name}/POSCAR"
        initSlab = read(initPoscar)

        fileName = "CONTCAR_" + name
        first_name = CONTCAR_DIRECTORY.split("/")[1]

        slab = read(f"{CONTCAR_DIRECTORY}/{fileName}")

        if not os.path.exists(f"images/{first_name}"):
            os.makedirs(f"images/{first_name}")
        if not os.path.exists(f"images/{POSCAR_DIRECTORY}_POSCAR"):
            os.makedirs(f"images/{POSCAR_DIRECTORY}_POSCAR")

        if os.path.exists(f"images/{first_name}/{name}"):
            if override:
                shutil.rmtree(f"images/{first_name}/{name}")
                shutil.rmtree(f"images/{POSCAR_DIRECTORY}_POSCAR/{name}")
            else:
                if os.path.exists(f"images/{first_name}/{name}/slab_135x_90y_225z.png"):
                    initImages1.append(
                        os.path.abspath(
                            f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_135x_90y_225z.png"
                        )
                    )
                    initImages2.append(
                        os.path.abspath(
                            f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_180x_180y_45z.png"
                        )
                    )
                    initImages3.append(
                        os.path.abspath(
                            f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_225x_225y_35z.png"
                        )
                    )
                    images1.append(
                        os.path.abspath(
                            f"images/{first_name}/{name}/slab_135x_90y_225z.png"
                        )
                    )
                    images2.append(
                        os.path.abspath(
                            f"images/{first_name}/{name}/slab_180x_180y_45z.png"
                        )
                    )
                    images3.append(
                        os.path.abspath(
                            f"images/{first_name}/{name}/slab_225x_225y_35z.png"
                        )
                    )
                    continue

        os.makedirs(f"images/{first_name}/{name}")
        os.makedirs(f"images/{POSCAR_DIRECTORY}_POSCAR/{name}")

        plotThenSaveAtoms(
            initSlab,
            135,
            90,
            225,
            ax,
            f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_135x_90y_225z.png",
        )
        plotThenSaveAtoms(
            initSlab,
            180,
            180,
            45,
            ax,
            f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_180x_180y_45z.png",
        )
        plotThenSaveAtoms(
            initSlab,
            225,
            225,
            35,
            ax,
            f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_225x_225y_35z.png",
        )

        plotThenSaveAtoms(
            slab, 135, 90, 225, ax, f"images/{first_name}/{name}/slab_135x_90y_225z.png"
        )
        plotThenSaveAtoms(
            slab, 180, 180, 45, ax, f"images/{first_name}/{name}/slab_180x_180y_45z.png"
        )
        plotThenSaveAtoms(
            slab, 225, 225, 35, ax, f"images/{first_name}/{name}/slab_225x_225y_35z.png"
        )

        initImages1.append(
            os.path.abspath(
                f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_135x_90y_225z.png"
            )
        )
        initImages2.append(
            os.path.abspath(
                f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_180x_180y_45z.png"
            )
        )
        initImages3.append(
            os.path.abspath(
                f"images/{POSCAR_DIRECTORY}_POSCAR/{name}/slab_225x_225y_35z.png"
            )
        )

        images1.append(
            os.path.abspath(f"images/{first_name}/{name}/slab_135x_90y_225z.png")
        )
        images2.append(
            os.path.abspath(f"images/{first_name}/{name}/slab_180x_180y_45z.png")
        )
        images3.append(
            os.path.abspath(f"images/{first_name}/{name}/slab_225x_225y_35z.png")
        )

    plt.close(fig)

    initImages1.sort()
    initImages2.sort()
    initImages3.sort()

    images1.sort()
    images2.sort()
    images3.sort()

    df["initialAngle1"] = initImages2
    df["initialAngle2"] = initImages1
    df["initialAngle3"] = initImages3

    df["finalAngle1"] = images2
    df["finalAngle2"] = images1
    df["finalAngle3"] = images3

    image_cols = [
        "initialAngle1",
        "initialAngle2",
        "initialAngle3",
        "finalAngle1",
        "finalAngle2",
        "finalAngle3",
    ]

    format_dict = {}
    for image_col in image_cols:
        format_dict[image_col] = path_to_image_html

    return df, format_dict


def getInitialXYfromDfAtoms(df, symbol: str, key: str, slab):
    xypairs = []
    atoms = getSurfaceAtoms(symbol, 0, slab)
    for name in df[key]:
        if len(name) == 2:
            index = name[1]
            x = atoms[int(index)].position[0]
            y = atoms[int(index)].position[1]
        else:
            indices = name.split(symbol)[1]
            x, y = find_average_of_symbol(
                symbol, [int(idx) for idx in list(indices)], slab, -2
            )

        xypairs.append((x, y))
    return xypairs


def addShortestThreeBondLengthsToDf(
    df, key: str, symbol1: str, symbol2: str, directory: str, starting: str
):
    names = df[key]
    formatted_list = []
    for name in names:
        fileName = starting + "_" + name
        slab = read(f"{directory}/{fileName}")
        _, dis = calculateDistancesForEachAtomPair(slab.copy(), symbol1, symbol2)
        formatted = f"{dis[0]}<br>{dis[1]}<br>{dis[2]}"
        formatted_list.append(formatted)

    refKey = f"Shortest distances between atoms of {symbol1}, {symbol2} (Å)"
    df[refKey] = formatted_list
    return refKey


slab = read("CNST_CONTCAR_WO3_T")
middle_slab = read("CNST_CONTCAR_WO3_M", format="vasp")
large_slab = read("CNST_CONTCAR_WO3_L", format="vasp")
backup_slab = read("backupPSCR", format="vasp")
emptyCell = read("CNST_CONTCAR_EMPTY")

old_height_above_slab = 2.2
height_above_slab = 1.5
height_above_slab_for_vacancies = 0.5
height_above_slab_for_vacancies_2 = 1
height_above_slab_for_H2_second_layer = -0.35

triangle_1 = [0, 1, 4]
triangle_2 = [2, 3, 5]

cleanUp()

from energies import h2_wo3_energy, wo3_v_energy, h2o_energy

# generateSlabVac(middle_slab.copy(), "O", 0)
# add_h2o_vacancy(
#     backup_slab.copy(), h2o.copy(), height_above_slab_for_vacancies, "O", 0, "O_down"
# )
# add_h2o_vacancy(
#     backup_slab.copy(), h2o.copy(), height_above_slab_for_vacancies, "O", 1, "O_down"
# )

# print(h2o_3_amount_energy - (h2o_energy + h2o_3_1vac_O0_energy))
# print(h2o_2_amount_energy - (h2o_energy + h2o_2_1vac_O0_energy))
# print(h2_wo3_energy - (h2o_energy + wo3_v_energy))
# print(h2_avg_o014_wo3_energy - (wo3_energy + h2_energy))
# # print(h2_bridge_wo3_energy - (wo3_energy + h2_energy))

# add_n2_vacancy(middle_slab.copy(), n2.copy(), height_above_slab_for_vacancies, "O", 0)

# print(medium_n2_vac_energy - (n2_energy + medium_wo3_v_energy))
# print(large_n2_vac_energy - (n2_energy + large_wo3_v_energy))
# print(large_n2_vac_energy)
# print(large_wo3_v_energy)

# print(medium_n2_vac_energy)
# print(medium_wo3_v_energy)

# print(n2_energy)

# import matplotlib.pyplot as plt
# from pymatgen.io.vasp import Vasprun

# sigma = 0.1  # smooth out the DOS
# vr = Vasprun("../001_WO3/vasprun.xml")

# dos = vr.complete_dos
# dos.densities = dos.get_smeared_densities(0.1)

# x = dos.energies - dos.efermi
# y = dos.get_densities()
# plt.plot(x, y, color="k")
# plt.xlabel("E - Ef (eV)")
# plt.ylabel("DOS (a.u.)")
# plt.show()


data = parseACFdat(slab)
print(data)

avgs = {}
for key in data.keys():
    item = data[key]
    avgs[key] = sum(item) / len(item)

print(avgs)


print("----done----")
