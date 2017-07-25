#!/bin/bash

# Get line number of delR
startLine=$(grep -n delR data | sed 's/:.*//g')
endLine=$(grep -n -m 1 -e '/' -e '=' data | sed 's:.*//g')

# Delete until backslash

