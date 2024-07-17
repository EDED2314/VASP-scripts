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

## Local Script Workspace Tree
```
.
├── Archive_WO3_b
│   ├── 2_H
│   │   └── O0
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── H
│   │   ├── Avg-O014
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── Avg-O235
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O0
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O1
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O2
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O3
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O4
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O5
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── W0
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── W1
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   └── W2
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── H2O
│   │   ├── V-O0-CD
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-CL
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-CR
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-CU
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-H2D
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-HDD
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-HDL
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-HDR
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-HDU
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   └── V-O0-OD
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── H2O-2H-WO3-V
│   │   └── O0
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── H_x2y2
│   │   └── O0
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── N2
│   │   ├── V-O0-CD
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-CL
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-CR
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── V-O0-CU
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   └── V-O0-UPR
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── N2_x2y2_V
│   │   ├── N_POTCAR
│   │   └── V-O0-UPR
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   ├── OSZICAR_WO3_V
│   ├── POST
│   │   ├── H2O_V
│   │   │   └── OSZICAR_EXTREMELY_CLOSE
│   │   ├── H_POST
│   │   │   ├── OSZICAR_Avg-O014
│   │   │   ├── OSZICAR_Avg-O235
│   │   │   ├── OSZICAR_O0
│   │   │   ├── OSZICAR_O1
│   │   │   ├── OSZICAR_O2
│   │   │   ├── OSZICAR_O3
│   │   │   ├── OSZICAR_O4
│   │   │   ├── OSZICAR_O5
│   │   │   ├── OSZICAR_W0
│   │   │   ├── OSZICAR_W1
│   │   │   └── OSZICAR_W2
│   │   └── N2_V
│   │       └── OSZICAR_N2_UPR_WO3_V
│   ├── POSTCONTCAR
│   │   └── H_CONTCAR
│   │       ├── CONTCAR_Avg-O014
│   │       ├── CONTCAR_Avg-O235
│   │       ├── CONTCAR_O0
│   │       ├── CONTCAR_O1
│   │       ├── CONTCAR_O2
│   │       ├── CONTCAR_O3
│   │       ├── CONTCAR_O4
│   │       ├── CONTCAR_O5
│   │       ├── CONTCAR_W0
│   │       ├── CONTCAR_W1
│   │       └── CONTCAR_W2
│   └── surface
│       ├── INCAR
│       ├── KPOINTS
│       ├── POSCAR
│       ├── POTCAR
│       └── gpu.slurm
├── CNST_CONTCAR_EMPTY
├── CNST_CONTCAR_WO3_B
├── CNST_CONTCAR_WO3_T
├── CNST_PSCR_WO3_LARGE
├── H
│   ├── 1stLayer
│   │   ├── O0
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O1
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   ├── O2
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   ├── POTCAR
│   │   │   └── gpu.slurm
│   │   └── P0.0
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       ├── POTCAR
│   │       └── gpu.slurm
│   └── 2ndLayer
│       ├── O0
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   ├── POTCAR
│       │   └── gpu.slurm
│       ├── O1
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   ├── POTCAR
│       │   └── gpu.slurm
│       ├── O2
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   ├── POTCAR
│       │   └── gpu.slurm
│       ├── O3
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   ├── POTCAR
│       │   └── gpu.slurm
│       ├── O4
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   ├── POTCAR
│       │   └── gpu.slurm
│       └── O5
│           ├── INCAR
│           ├── KPOINTS
│           ├── POSCAR
│           ├── POTCAR
│           └── gpu.slurm
├── H2O
│   ├── V-O0-OD
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   ├── V-O1-OD
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   ├── V-O2-OD
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   └── V-P0.0-OD
│       ├── INCAR
│       ├── KPOINTS
│       ├── POSCAR
│       ├── POTCAR
│       └── gpu.slurm
├── KPOINTS_POSCAR
├── N2
│   ├── N_POTCAR
│   ├── V-O0-UPR
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   ├── V-O1-UPR
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   └── V-O2-UPR
│       ├── INCAR
│       ├── KPOINTS
│       ├── POSCAR
│       ├── POTCAR
│       └── gpu.slurm
├── OSZICAR_PLOTTER.py
├── POSTOUTPUT
│   ├── OSZICAR_H
│   ├── OSZICAR_H2
│   ├── OSZICAR_H2O
│   ├── OSZICAR_N
│   ├── OSZICAR_N2
│   └── OSZICAR_WO3
├── ads.py
├── ads.zip
├── adsorbates
│   ├── H
│   │   ├── INCAR_H
│   │   ├── KPOINTS
│   │   ├── POSCAR_H
│   │   ├── POTCAR_H
│   │   └── gpu.slurm
│   ├── H2
│   │   ├── INCAR_H2
│   │   ├── KPOINTS
│   │   ├── POSCAR_H2
│   │   ├── POTCAR_H2
│   │   └── gpu.slurm
│   ├── H2O
│   │   ├── INCAR_H2O
│   │   ├── KPOINTS
│   │   ├── POSCAR_H2O
│   │   ├── POTCAR_H2O
│   │   └── gpu.slurm
│   └── N2
│       ├── INCAR_N2
│       ├── KPOINTS
│       ├── POSCAR_N2
│       ├── POTCAR_N2
│       └── gpu.slurm
├── backupPSCR
├── data
│   ├── H_atom_adsorption_energy.html
│   ├── adsorption_energy.csv
│   ├── slab_135x_90y_225z.png
│   ├── slab_180x_180y_45z.png
│   └── slab_225x_225y_35z.png
├── images
│   ├── H_CONTCAR
│   │   ├── Avg-O014
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── Avg-O235
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── O0
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── O1
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── O2
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── O3
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── O4
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── O5
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── W0
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   ├── W1
│   │   │   ├── slab_135x_90y_225z.png
│   │   │   ├── slab_180x_180y_45z.png
│   │   │   └── slab_225x_225y_35z.png
│   │   └── W2
│   │       ├── slab_135x_90y_225z.png
│   │       ├── slab_180x_180y_45z.png
│   │       └── slab_225x_225y_35z.png
│   └── H_POSCAR
│       ├── Avg-O014
│       │   ├── slab_135x_90y_225z.png
│       │   ├── slab_180x_180y_45z.png
│       │   └── slab_225x_225y_35z.png
│       ├── O0
│       │   ├── slab_135x_90y_225z.png
│       │   ├── slab_180x_180y_45z.png
│       │   └── slab_225x_225y_35z.png
│       ├── W0
│       │   ├── slab_135x_90y_225z.png
│       │   ├── slab_180x_180y_45z.png
│       │   └── slab_225x_225y_35z.png
│       └── W1
│           ├── slab_135x_90y_225z.png
│           ├── slab_180x_180y_45z.png
│           └── slab_225x_225y_35z.png
├── notebooks
│   ├── gpaw.txt
│   ├── helpme.traj
│   └── main.ipynb
├── requirements.txt
├── surface
│   ├── INCAR
│   ├── KPOINTS
│   ├── POSCAR
│   ├── POTCAR
│   └── gpu.slurm
├── surface_v
│   ├── O0
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   ├── O1
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   ├── POTCAR
│   │   └── gpu.slurm
│   └── O2
│       ├── INCAR
│       ├── KPOINTS
│       ├── POSCAR
│       ├── POTCAR
│       └── gpu.slurm
├── surface_v.zip
├── surface_x2y2
│   ├── INCAR
│   ├── KPOINTS
│   ├── POSCAR
│   ├── POTCAR
│   └── gpu.slurm
├── templates_W001
│   ├── INCAR
│   ├── KPOINTS
│   ├── POTCAR
│   └── gpu.slurm
├── templates_W001_x2y2
│   ├── INCAR
│   ├── KPOINTS
│   ├── POTCAR
│   └── gpu.slurm
├── templates_adsorbate
│   ├── INCAR_H
│   ├── INCAR_H2
│   ├── INCAR_H2O
│   ├── INCAR_N2
│   ├── KPOINTS
│   ├── POTCAR_H
│   ├── POTCAR_H2
│   ├── POTCAR_H2O
│   ├── POTCAR_N2
│   ├── POTCAR_O
│   └── gpu.slurm
└── vis.ipynb
```

## Super Computer Reference Workspace Tree
```
.
├── ads
│   ├── H
│   │   ├── 1stLayer
│   │   │   ├── O0
│   │   │   │   ├── gpu.slurm
│   │   │   │   ├── INCAR
│   │   │   │   ├── KPOINTS
│   │   │   │   ├── POSCAR
│   │   │   │   └── POTCAR
│   │   │   ├── O1
│   │   │   │   ├── gpu.slurm
│   │   │   │   ├── INCAR
│   │   │   │   ├── KPOINTS
│   │   │   │   ├── POSCAR
│   │   │   │   └── POTCAR
│   │   │   ├── O2
│   │   │   │   ├── gpu.slurm
│   │   │   │   ├── INCAR
│   │   │   │   ├── KPOINTS
│   │   │   │   ├── POSCAR
│   │   │   │   └── POTCAR
│   │   │   └── P0.0
│   │   │       ├── gpu.slurm
│   │   │       ├── INCAR
│   │   │       ├── KPOINTS
│   │   │       ├── POSCAR
│   │   │       └── POTCAR
│   │   └── 2ndLayer
│   │       ├── O0
│   │       │   ├── gpu.slurm
│   │       │   ├── INCAR
│   │       │   ├── KPOINTS
│   │       │   ├── POSCAR
│   │       │   └── POTCAR
│   │       ├── O1
│   │       │   ├── gpu.slurm
│   │       │   ├── INCAR
│   │       │   ├── KPOINTS
│   │       │   ├── POSCAR
│   │       │   └── POTCAR
│   │       ├── O2
│   │       │   ├── gpu.slurm
│   │       │   ├── INCAR
│   │       │   ├── KPOINTS
│   │       │   ├── POSCAR
│   │       │   └── POTCAR
│   │       ├── O3
│   │       │   ├── gpu.slurm
│   │       │   ├── INCAR
│   │       │   ├── KPOINTS
│   │       │   ├── POSCAR
│   │       │   └── POTCAR
│   │       ├── O4
│   │       │   ├── gpu.slurm
│   │       │   ├── INCAR
│   │       │   ├── KPOINTS
│   │       │   ├── POSCAR
│   │       │   └── POTCAR
│   │       └── O5
│   │           ├── gpu.slurm
│   │           ├── INCAR
│   │           ├── KPOINTS
│   │           ├── POSCAR
│   │           └── POTCAR
│   ├── H2O
│   │   ├── V-O0-OD
│   │   │   ├── gpu.slurm
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   └── POTCAR
│   │   ├── V-O1-OD
│   │   │   ├── gpu.slurm
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   └── POTCAR
│   │   ├── V-O2-OD
│   │   │   ├── gpu.slurm
│   │   │   ├── INCAR
│   │   │   ├── KPOINTS
│   │   │   ├── POSCAR
│   │   │   └── POTCAR
│   │   └── V-P0.0-OD
│   │       ├── gpu.slurm
│   │       ├── INCAR
│   │       ├── KPOINTS
│   │       ├── POSCAR
│   │       └── POTCAR
│   └── N2
│       ├── N_POTCAR
│       ├── V-O0-UPR
│       │   ├── gpu.slurm
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   └── POTCAR
│       ├── V-O1-UPR
│       │   ├── gpu.slurm
│       │   ├── INCAR
│       │   ├── KPOINTS
│       │   ├── POSCAR
│       │   └── POTCAR
│       └── V-O2-UPR
│           ├── gpu.slurm
│           ├── INCAR
│           ├── KPOINTS
│           ├── POSCAR
│           └── POTCAR
├── surface_v
│   ├── O0
│   │   ├── gpu.slurm
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   └── POTCAR
│   ├── O1
│   │   ├── gpu.slurm
│   │   ├── INCAR
│   │   ├── KPOINTS
│   │   ├── POSCAR
│   │   └── POTCAR
│   └── O2
│       ├── gpu.slurm
│       ├── INCAR
│       ├── KPOINTS
│       ├── POSCAR
│       └── POTCAR
└── vac
```
