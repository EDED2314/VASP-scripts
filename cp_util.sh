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

# Prompt the user for the file type
read -p "Enter the file type to search for (e.g., XDATCAR, OSZICAR, etc.): " FILE_TYPE

# Prompt the user for the directory
read -p "Enter the directory name (e.g., H): " SOURCE_DIR

# Change this for the folder that saves the oszicar files
DEST_DIR="${SOURCE_DIR}_${FILE_TYPE}"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all files of the specified type in the source directory and its subdirectories
find "$SOURCE_DIR" -type f -name "$FILE_TYPE" | while read -r FILE; do
    # Extract the subdirectory path relative to the source directory
    RELATIVE_PATH=$(dirname "${FILE#$SOURCE_DIR/}")
    # Create the corresponding subdirectory in the destination directory
    mkdir -p "$DEST_DIR/$RELATIVE_PATH"
    # Copy the file to the destination directory, preserving the relative path
    cp "$FILE" "$DEST_DIR/$RELATIVE_PATH"
    echo "Copied $FILE to $DEST_DIR/$RELATIVE_PATH"
done

# Zip the destination directory and remove the original directory
zip -r "${DEST_DIR}.zip" "$DEST_DIR"
rm -rf "$DEST_DIR"

echo "done"
echo "Zip file location: ${PWD}/${DEST_DIR}.zip"



