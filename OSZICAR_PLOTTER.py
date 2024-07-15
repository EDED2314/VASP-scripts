#!/usr/bin/python
# -*- coding:utf-8 -*-
# Re-purposed from https://github.com/K4ys4r/VASP_Scripts/blob/master/OSZICAR_MD_Plot/Plot_OSZICAR_MD.py

import pylab as p

ff = " as a function of MD Steps"

"""
The next line (after the N+1 lines) gives information about the total energy 
after obtaining convergence. The first values is the total free energy F 
(at this point the energy of the reference atom has been subtracted),
E0 is the energy for sigma → 0 (see also Partial occupancies) 
and dE is the change in the total energy between the current and the last step; 
for a static run dE is the entropy multiplied by  sigma.
"""

OSZ_Labels = {
    "Steps": [0, "Step ", " (time) "],
    "F": [1, "Total Free Energy", " (eV) ", "Free Energy" + ff],
    "E0": [2, "Energy sigma -> 0", " (eV) ", "Energy" + ff],
    "dE": [3, "Change in total energy", " (eV) ", "Free Energy" + ff],
}


# Molecular Dynamic OSZICAR Reading
def OSZICAR_READ(fileName):
    inp = open("POSTOUTPUT/" + fileName, "r")
    f = inp.readlines()
    inp.close()
    DATA = p.array([])
    for i in range(len(f)):
        if "F=" in f[i]:
            info = f[i].split()
            step = int(info[0])
            if step < 5:
                # this causes some errors for sims that have less than 5 steps
                continue
            info.pop(5)
            info[5] = "dE="
            info[6] = info[6].strip("=")
            # print(info)
            if len(DATA) == 0:
                DATA = p.array([info[::2]], dtype=float)
            else:
                DATA = p.append(DATA, [p.array(info[::2], dtype=float)], axis=0)
    return DATA


##


# Plot Data
def PLOT_DATA(arr, Xplot, Yplot):
    x = arr[:, OSZ_Labels[Xplot][0]]
    y = arr[:, OSZ_Labels[Yplot][0]]
    lb = OSZ_Labels[Yplot][1]
    print(
        "\n\tMean {} = {} {}\n".format(
            OSZ_Labels[Yplot][1], p.mean(y), OSZ_Labels[Yplot][2]
        )
    )
    p.rcParams["font.family"] = "serif"
    p.figure(figsize=(10, 6))
    p.rcParams.update({"font.size": 14})
    p.xlabel(OSZ_Labels[Xplot][1] + OSZ_Labels[Xplot][2])
    p.ylabel(OSZ_Labels[Yplot][1] + OSZ_Labels[Yplot][2])
    p.title(OSZ_Labels[Yplot][3])
    p.plot(x, y, "r")
    p.show()


ending = input("File ending\n>>>")
OSZ = OSZICAR_READ("OSZICAR_" + ending)
# print(OSZ)

Xplot = "Steps"
Yplot = input("Yplot (F, E0, dE ) > ")

PLOT_DATA(OSZ, Xplot, Yplot)
