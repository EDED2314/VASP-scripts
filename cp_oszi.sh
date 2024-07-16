#!/bin/bash

# TO RUN
# bash cp_oszi.sh

#Run in the directory that has a file structure as such:
#  | Current Directory
#  |
#  | Subfolder 1
#  | -----  | VASP OUTPUT FILES....
#  | -----  | POTCAR
#  | -----  | OSZICAR
#  | -----  | POSCAR ....
#  | -----  | VASP OUTPUT FILES....
#  | Subfolder 2
#  | -----  | VASP OUTPUT FILES.......
#  | Subfolder 3
#  | -----  | VASP OUTPUT FILES..........
#  | Subfolder X....


# Define the source and destination directories
SOURCE_DIR="H"

#Change this for the folder that saves the oszicar files
DEST_DIR="H_POST"


# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through each subdirectory in the source directory
for SUBDIR in "$SOURCE_DIR"/*; do
    if [ -d "$SUBDIR" ]; then
        # Check if the OSZICAR file exists in the subdirectory
        if [ -f "$SUBDIR/OSZICAR" ]; then
            # Extract the subdirectory name
            SUBDIR_NAME=$(basename "$SUBDIR")
            # Copy and rename the OSZICAR file to the destination directory
            cp "$SUBDIR/OSZICAR" "$DEST_DIR/OSZICAR_${SUBDIR_NAME}"
            echo "Copied and renamed $SUBDIR/OSZICAR to $DEST_DIR/OSZICAR_${SUBDIR_NAME}"
        else
            echo "OSZICAR not found in $SUBDIR"
        fi
    fi
done


zip -r "$DEST_DIR".zip $DEST_DIR
rm -rf $DEST_DIR

