#!/bin/sh
# Check if an argument is provided
if [ $# -eq 0 ]; then
    echo "Usage: $0 <your_argument1> <your_argument2>"
    exit 1
fi

# Access the first argument using $1
your_argument="$1"
your_argument2="$2"

if [ -z "$2" ]; then
    your_argument2="FALSE"  # Set your default value here
else
    your_argument2="$2"
fi

echo "Sample Name: $your_argument"
echo "Combos: $your_argument2"

export MY_ARGUMENT="$your_argument"
export MY_ARGUMENT2="$your_argument2"

Rscript Tools/makeDB.R 
