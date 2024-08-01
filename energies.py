from constants import *


def parseACFdat(slab):
    lines = []
    with open("acf.dat", "r") as f:
        lines = f.readlines()

    cleaned_lines = []
    for line in lines:
        if not ("------" in line or "#" in line or "V" in line or "NUM" in line):
            cleaned_lines.append(" ".join(line.strip("\n").split()).split(" "))

    atomIdx = []
    charge = []
    for line in cleaned_lines:
        atomIdx.append(line[0])
        charge.append(line[4])

    symbols = []
    z_vals = []
    for atom in slab:
        sym = atom.symbol
        if sym == "W":
            z_vals.append(W)
        elif sym == "O":
            z_vals.append(O)
        elif sym == "H":
            z_vals.append(H)
        elif sym == "N":
            z_vals.append(N)

        symbols.append(sym)

    assert len(z_vals) == len(symbols) == len(atomIdx) == len(charge)

    atomicCharge = np.subtract(
        np.array(charge, dtype="float64"), np.array(z_vals, dtype="float64")
    )

    finalData = {}
    for i in range(len(list(atomicCharge))):
        if not symbols[i] in finalData:
            finalData[symbols[i]] = []
        finalData[symbols[i]].append(atomicCharge[i])

    return finalData


def readOszicarFileAndGetLastLineEnergy(fileName: str, debug=False):
    lines = []
    with open(fileName) as f:
        lines = f.readlines()

    lastLine = lines[-1]
    tmp = lastLine.split()
    energy = float(tmp[2])

    return energy


W = 14.000
O = 6.000
H = 1.000
N = 5.000

o_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_O")
o2_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_O2")
no_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_NO")
n2o_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_N2O")
n_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_N")
n2_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_N2")
h2o_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_H2O")
h2_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_H2")
h_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_H")
wo3_energy = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_WO3")
wo3_v_energy = (
    readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_WO3_V_O0")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_WO3_V_O1")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/OSZICAR_WO3_V_O2")
) / 3.0

n2_vac_energy = (
    readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N2_OSZICAR/OSZICAR_V-O0-UPR")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N2_OSZICAR/OSZICAR_V-O1-UPR")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N2_OSZICAR/OSZICAR_V-O2-UPR")
) / 3.0

n_energy_O0 = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N_OSZICAR/O0_NC")
n2_energy_O0 = readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N2_OSZICAR/O0")

n_vac_energy = (
    readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N_OSZICAR/OSZICAR_O2")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/N_OSZICAR/OSZICAR_O0")
) / 2.0
h_wo3_energy = (
    readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/H_1stLayer_OSZICAR/OSZICAR_O0")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/H_1stLayer_OSZICAR/OSZICAR_O1")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/H_1stLayer_OSZICAR/OSZICAR_O2")
) / 3.0
h2_wo3_energy = (
    readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/H2O_OSZICAR/OSZICAR_V-O0-OD")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/H2O_OSZICAR/OSZICAR_V-O1-OD")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/H2O_OSZICAR/OSZICAR_V-O2-OD")
) / 3.0

h2_avg_o014_wo3_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/H2_OSZICAR/avgO014"
)
h2_bridge_wo3_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/H2_OSZICAR/bridge0"
)

h2o_2_amount_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/H2O_amt_OSZICAR/OSZICAR_2H2O"
)
h2o_3_amount_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/H2O_amt_OSZICAR/OSZICAR_3H2O"
)
h2o_2_1vac_O0_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/H2O_amt_OSZICAR/OSZICAR_1VAC_2H2O"
)
h2o_3_1vac_O0_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/H2O_amt_OSZICAR/OSZICAR_1VAC_3H2O"
)


large_wo3_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/Large/OSZICAR_WO3"
)
large_wo3_v_energy = (
    readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/Large/OSZICAR_WO3_V_O1")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/Large/OSZICAR_WO3_V_O2")
    + readOszicarFileAndGetLastLineEnergy(f"{OUTPUT_DIR}/Large/OSZICAR_WO3_V_O3")
) / 3.0

large_n2_vac_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/Large/N2_OSZICAR/O1"
)


medium_wo3_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/Medium/OSZICAR_WO3"
)
medium_wo3_v_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/Medium/OSZICAR_WO3_V_O0"
)
medium_n2_vac_energy = readOszicarFileAndGetLastLineEnergy(
    f"{OUTPUT_DIR}/Medium/N2_OSZICAR/O0_VAC"
)
