#!/bin/bash

e3post.py --job=binary_diffusion  --slice-list="0,:,5,0" --output-file=numerical-profile-Ficks.data

echo "----------------------------------------------------------"
echo "Check simulation against exact solution:"
e3post.py --job=binary_diffusion --ref-function=ref_function_wrapper.py
