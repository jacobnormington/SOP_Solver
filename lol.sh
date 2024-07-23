#!/bin/bash

# Define the target value to compare against
TARGET=71749

while true; do
    # Run the command and capture the output
    output=$(./sop_solver soplib/R.200.100.60.sop 32 soplib_config.txt)
    
    # Extract the last line of the output
    last_line=$(echo "$output" | tail -n 1)
    
    # Extract the first number from the last line (format: x,y)
    result=$(echo "$last_line" | cut -d ',' -f 1)
    
    # Print the output for debugging purposes (optional)
   # echo "Output: $output"
    echo "Last Line: $last_line"
   # echo "Result: $result"
    
    # Check if the result is not equal to the target value
    if [ "$result" -ne "$TARGET" ]; then
        echo "Found a result not equal to $TARGET: $result"
        break
    fi
done