#/usr/bin/env bash

if [ -z "$1" ]; then
  echo "Parses a directory with text files containg /usr/bin/time outputs\nUsage: $0 <osx|gnu> <folder>"
  exit 1
fi


# Check if a folder argument is provided
if [ -z "$2" ]; then
  echo "Parses a directory with text files containg /usr/bin/time outputs\nUsage: $0 <osx|gnu> <folder>"
  exit 1
fi

CHOICE="$1"

# Assign the folder to a variable
FOLDER="$2"

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
  name=$(echo "$filename" | sed -E 's/^(.*)\.bm\..*\.time\.txt$/\1/')
  case "$CHOICE" in

    osx)
      # Use grep to extract lines containing the keyword "time"
      timing=$(grep -i "real" "$file" | sed -E 's/^[[:space:]]*([0-9]+,[0-9]+)[[:space:]]+real.*/\1/')
      timing=$(echo "$timing" | sed -E 's/,/./g')
      memory=$(grep -i "peak memory" "$file" | sed -E 's/^[[:space:]]*([0-9]+)[[:space:]]+peak.*/\1/') 
      echo "${name}\t${tool_name}\t${timing}\t${memory}"
      ;;

    gnu)# Extract timing using sed
      timing=$(grep -i "Elapsed" "$file" | sed -E 's/.*time \(.*\): ([0-9]+:[0-9]+(\.[0-9]+|:[0-9]+)).*/\1/')

      # Check if the extracted timing is in mm:ss.ss or hh:mm:ss format
      if [[ "$timing" == *.* ]]; then
        # Case: mm:ss.ss format
        minutes=$(echo "$timing" | sed -E 's/:.*//') # Extract minutes
        seconds=$(echo "$timing" | sed -E 's/^[0-9]+://') # Extract seconds
        total_seconds=$(echo "$minutes * 60 + $seconds" | bc) # Calculate total seconds
      else
        # Case: hh:mm:ss format
        hours=$(echo "$timing" | sed -E 's/:.*//') # Extract hours
        minutes=$(echo "$timing" | sed -E 's/^[0-9]+:([0-9]+):.*/\1/') # Extract minutes
        seconds=$(echo "$timing" | sed -E 's/^.*://') # Extract seconds
        total_seconds=$(echo "$hours * 3600 + $minutes * 60 + $seconds" | bc) # Calculate total seconds
      fi


      memory=$(grep -i "Maximum resident" "$file" | sed -E 's/.* ([0-9]+)$/\1/') 
      echo "${name}\t${tool_name}\t${total_seconds}\t${memory}"
      ;;

    *)
      echo "Invalid option: $choice"
      echo "Usage: $0 <osx|gnu>"
      exit 1
      ;;
  esac

done
