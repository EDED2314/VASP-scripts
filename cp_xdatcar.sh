#!/bin/bash

# Define the source and destination directories
SOURCE_DIR="H"
DEST_DIR="H_XDATCAR"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Loop through each subdirectory in the source directory
for SUBDIR in "$SOURCE_DIR"/*; do
    if [ -d "$SUBDIR" ]; then
        # Check if the OSZICAR file exists in the subdirectory
        if [ -f "$SUBDIR/XDATCAR" ]; then
            # Extract the subdirectory name
            SUBDIR_NAME=$(basename "$SUBDIR")
            # Copy and rename the OSZICAR file to the destination directory
            cp "$SUBDIR/XDATCAR" "$DEST_DIR/XDATCAR_${SUBDIR_NAME}"
            echo "Copied and renamed $SUBDIR/XDATCAR to $DEST_DIR/XDATCAR_${SUBDIR_NAME}"
        else
            echo "XDATCAR not found in $SUBDIR"
        fi
    fi
done

zip -r "$DEST_DIR".zip $DEST_DIR
