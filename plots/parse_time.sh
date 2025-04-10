#/usr/bin/env bash

# Check if a folder argument is provided
if [ -z "$1" ]; then
    echo "Parses a directory with text files containg /usr/bin/time outputs\nUsage: $0 <folder>"
    exit 1
fi

# Assign the folder to a variable
FOLDER="$1"

# Check if the provided argument is a valid directory
if [ ! -d "$FOLDER" ]; then
    echo "Error: $FOLDER is not a valid directory."
    exit 1
fi

# Process files ending with *.time.txt in the folder
for file in "$FOLDER"/*.time.txt; do
    # Check if there are matching files (to avoid wildcard issues)
    if [ ! -e "$file" ]; then
        echo "No files ending in *.time.txt found in $FOLDER."
        exit 0
    fi
    filename=$(basename "$file") # Get the filename without the path
    tool_name=$(echo "$filename" | sed -E 's/.*\.([^.]+)\.time\.txt/\1/') # Extract tool_name using regexZZ


    
    # Use grep to extract lines containing the keyword "time"
    timing=$(grep -i "real" "$file" | sed -E 's/^[[:space:]]*([0-9]+,[0-9]+)[[:space:]]+real.*/\1/')
    timing=$(echo "$timing" | sed -E 's/,/./g')
    memory=$(grep -i "peak memory" "$file" | sed -E 's/^[[:space:]]*([0-9]+)[[:space:]]+peak.*/\1/') 
    name=$(echo "$filename" | sed -E 's/^(.*)\.bm\..*\.time\.txt$/\1/')
    echo "${name}\t${tool_name}\t${timing}\t${memory}"
done
