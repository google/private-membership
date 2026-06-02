#!/bin/bash

# A script to generate a large file where each line is unique.

# --- Configuration ---
# Number of lines (2^16)
LINE_COUNT=32768

# Number of characters per line
CHAR_COUNT=8192

# Output filename
OUTPUT_FILE="db-large.txt"

# --- Generation ---

echo "⏳ Generating $LINE_COUNT lines of $CHAR_COUNT characters each..."

# Ensure the output file is empty before starting
> "$OUTPUT_FILE"

# Loop LINE_COUNT times
for (( i=1; i<=$LINE_COUNT; i++ ))
do
  # 1. Generate one line of random alphanumeric characters.
  #    - 'tr' filters random bytes to keep only 'A-Za-z0-9'.
  #    - 'head -c' takes exactly CHAR_COUNT characters for the line.
  # 2. Append the result to the output file.
  #    - The command is run inside a command substitution '$(...)'
  #    - 'echo' prints the result and automatically adds the required newline.
  #    - '>>' appends the output to the file.
  echo "$(LC_ALL=C tr -dc 'A-Za-z0-9' < /dev/urandom | head -c $CHAR_COUNT)" >> "$OUTPUT_FILE"
done

echo "✅ Done. File '$OUTPUT_FILE' has been created."