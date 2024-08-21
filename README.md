# VASP (Vienna Ab initio Simulation Package) Scripts
More about VASP [here](https://www.vasp.at/)

## Features
- Initial structure gen for **POSCAR** and KPOINTS (POTCAR, INCAR, KPOINTS, and job file (.slurm) will need to be provided in the templates folder) 
- OSZICAR plotter (Y v.s. steps)
- Adsorption energy calculator
- Bond length calculations between atom pairs
- HTML webpage generation with data and important previews
- More to come

### Functions
- Freezing bottom two layers of atoms (can be modified to freeze bottom *n* layers)
- Simulation folder generation from the template folder
- Slabs, molecules, and atom manipulation (Adsorbates on a slab, Atom/Molecule in a vacuum, Vacancy on a slab)
- More to come

## HTML output example
![image](https://github.com/EDED2314/VASP-scripts/blob/main/HTML%20Output%20Example%207.16.24.jpg)

## Graph output example

### Coming soon (feature already done... don't know if I can upload the actual graph though)

## Local Workspace Reference Tree
```
.

├── H
│   ├── 1stLayer
│   │   ├── O0
│   │   ├── O1
│   │   ├── O2
│   │   └── P0.0
│   └── 2ndLayer
│       ├── O0
│       ├── O1
│       ├── O2
│       ├── O3
│       ├── O4
│       └── O5
├── H2O
│   ├── V-O0-OD
│   ├── V-O1-OD
│   └── V-O2-OD
├── H2O_amount
│   ├── 2-O0-O1
│   └── 3-O0-O1-O2
├── H2O_x1y2
│   └── V-O0-OD
├── H2O_x2y2
│   ├── V-O1-OD
│   ├── V-O2-OD
│   └── V-O3-OD
├── H_x1y2
│   ├── 1stL-O0
│   └── 2ndL-O0
├── H_x2y2
│   ├── O1
│   ├── O2
│   └── O3
├── N
│   ├── O0
│   └── O2
├── N2
│   ├── V-O0-UPR
│   ├── V-O1-UPR
│   └── V-O2-UPR
├── N2_x1y2
│   └── V-O0-UPR
├── N2_x2y2
│   ├── V-O0-UPR
│   ├── V-O1-UPR
│   └── V-O2-UPR
├── POSTCONTCAR
│   ├── H2O_CONTCAR
│   ├── H_1stLayer_CONTCAR
│   ├── H_2ndLayer_CONTCAR
│   └── N2_CONTCAR
├── POSTOUTPUT
│   ├── Large
│   │   ├── H
│   │   ├── H2O_VAC_OSZICAR
│   │   └── N2_VAC_OSZICAR
│   ├── Medium
│   │   ├── H
│   │   ├── H2O
│   │   ├── H2O_amt_OSZICAR
│   │   └── N2_OSZICAR
│   ├── Small
│   │   ├── H2O_OSZICAR
│   │   ├── H2O_amt_OSZICAR
│   │   ├── H2_OSZICAR
│   │   ├── H_1stLayer_OSZICAR
│   │   ├── N2_OSZICAR
│   │   └── N_OSZICAR
│   ├── gas
│   ├── spin_Large
│   ├── spin_Medium
│   └── spin_Small
│       ├── H
│       │   ├── 1LO0
│       │   └── 2LO0
│       ├── H2O
│       │   └── V-O0-OD
│       ├── H2O_amount
│       │   ├── 2-O0-O1
│       │   ├── 2-vac-O0
│       │   ├── 3-O0-O1-O2
│       │   └── 3-vac-O0
│       ├── N
│       │   ├── O0
│       │   └── vO0
│       └── N2
│           ├── O0
│           ├── V-O0-UPR
│           └── bridge0
├── ads_x1y2
│   ├── H2O_x1y2
│   │   └── V-O0-OD
│   ├── H_x1y2
│   │   ├── 1stL-O0
│   │   └── 2ndL-O0
│   └── N2_x1y2
│       └── V-O0-UPR
├── ads_x2y2
│   ├── H2O_x2y2
│   │   ├── V-O1-OD
│   │   ├── V-O2-OD
│   │   └── V-O3-OD
│   ├── H_x2y2
│   │   ├── O1
│   │   ├── O2
│   │   ├── O3
│   │   ├── V-O1-OD
│   │   ├── V-O2-OD
│   │   └── V-O3-OD
│   └── N2_x2y2
│       ├── V-O0-UPR
│       ├── V-O1-UPR
│       └── V-O2-UPR
├── adsorbates
│   ├── H
│   ├── H2
│   ├── H2O
│   ├── N
│   ├── N2
│   ├── N2O
│   ├── NO
│   └── O2
├── notebooks
├── surface
├── surface_v
│   ├── O0
│   ├── O1
│   └── O2
├── surface_x1y2
├── surface_x2y2
├── templates_W001
├── templates_W001_x1y2
├── templates_W001_x2y2
└── templates_adsorbate
```
## Super Computer Workspace Reference Tree
```
.
├── ads
│   ├── H
│   │   ├── 1stLayer
│   │   │   ├── O0
│   │   │   ├── O1
│   │   │   ├── O2
│   │   │   └── P0.0
│   │   └── 2ndLayer
│   │       ├── O0
│   │       ├── O1
│   │       ├── O2
│   │       ├── O3
│   │       ├── O4
│   │       └── O5
│   ├── H2O
│   │   ├── V-O0-OD
│   │   ├── V-O1-OD
│   │   ├── V-O2-OD
│   │   └── V-P0.0-OD
│   └── N2
│       ├── V-O0-UPR
│       ├── V-O1-UPR
│       └── V-O2-UPR
├── surface_v
│   ├── O0
│   ├── O1
│   └── O2
└── vac
```
