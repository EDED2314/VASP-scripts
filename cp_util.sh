#!/bin/bash

# TO RUN
# bash cp_util.sh

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

#!/bin/bash

# Prompt the user for the file type
read -p "Enter the file type to search for (e.g., XDATCAR, OSZICAR, etc.): " FILE_TYPE

#Prompt the user for the directory
read -p "Enter the directory name (e.g. H): " SOURCE_DIR

#Change this for the folder that saves the oszicar files
DEST_DIR="H_${FILE_TYPE}"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through each subdirectory in the source directory
for SUBDIR in "$SOURCE_DIR"/*; do
    if [ -d "$SUBDIR" ]; then
        # Check if the specified file type exists in the subdirectory
        if [ -f "$SUBDIR/$FILE_TYPE" ]; then
            # Extract the subdirectory name
            SUBDIR_NAME=$(basename "$SUBDIR")
            # Copy and rename the specified file to the destination directory
            cp "$SUBDIR/$FILE_TYPE" "$DEST_DIR/${FILE_TYPE}_${SUBDIR_NAME}"
            echo "Copied and renamed $SUBDIR/$FILE_TYPE to $DEST_DIR/${FILE_TYPE}_${SUBDIR_NAME}"
        else
            echo "$FILE_TYPE not found in $SUBDIR"
        fi
    fi
done

# Zip the destination directory and remove the original directory
zip -r "${DEST_DIR}.zip" "$DEST_DIR"
rm -rf "$DEST_DIR"

echo "done"
echo "Zip file location: ${PWD}/${DEST_DIR}.zip" 




