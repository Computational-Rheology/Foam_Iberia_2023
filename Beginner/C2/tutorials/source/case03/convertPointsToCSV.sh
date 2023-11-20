#!/bin/bash 

# Note:   Run this script from the main level of a case dir, or place the proper location
#     | of the blockMeshDict in the input_file var
#     |   The points.csv file must be place in the main level of a case dir or adjust the
#     | a01blockMeshView.pvsm file to read this file properly.


# User customization for the Input and output file names

  input_file="./system/blockMeshDict"
  output_file="points.csv"


# Script Section
# Create the header for the .csv file
echo "x,y,z,id" > "$output_file"

# Initialize the ID counter for the vertices' IDs
id=0
# Use a flag to identify the "vertices" section in the input file
inside_vertices=0


# Use grep to exclude lines that match the script itself
grep -v "^\s*# Loop through each line in the input file" "$input_file" | 
while IFS= read -r line; do
  # Check if the current line  in the while loop contains "vertices"
  if [[ $line == *vertices* ]]; then
    # Marking the start of the section when found
    inside_vertices=1
    continue
  fi
  # If we're inside the "vertices" section then check the line
  if [ "$inside_vertices" -eq 1 ]; then
    # Remove C++ style comments (//) and trim leading/trailing spaces
    cleaned_line=$(echo "$line" | sed 's|//.*||' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//')
    # Remove nested C|C++-style comments (/* ... */)
    cleaned_line=$(echo "$cleaned_line" | sed -E ':a; s|/\*[^*]*(\*[^/][^*]*)*\*/||g; ta')
    # Check if the cleaned line is not empty
    if [ -n "$cleaned_line" ]; then
      # Remove parentheses "(" and ")" from cleaned_line
      cleaned_line=$(echo "$cleaned_line" | tr -d '()')
      # Check if the cleaned line is not empty and not equal to anything with ;
      if [ -n "$cleaned_line" ] && [[ "$cleaned_line" != *";" ]]; then
        # Extract and print the coordinates with the current ID
        echo "$cleaned_line $id" | awk '{print $1","$2","$3","$4}' >> "$output_file"
        # Increment the ID
        ((id++))
      fi
    fi
    # Check if we've reached the end of the "vertices" section
    if [[ $line == *";" ]]; then
      inside_vertices=0
    fi
  fi
done

# Notify completion
echo "CSV file generated: $output_file"

# End WdCG