#!/bin/bash
### This scripts iterates BEAM_error_difference.py with different error amplitudes, given phase or amp option.
### Usage: bash run.sh phase 
###        bash run.sh amp

error="${1}" # phase or amplitude

array=(1.0 5.0 10.0 30.0 50.0 100.0)
for element in "${array[@]}"
do
    echo "amplitude = $element"
    if [ "$error" == "amp" ] #amplitude error
    then
        echo "amp"
        cmd="python BEAM_error_difference.py 8192 2048 10.0 1.0 3.0 y amp $element"
        echo $cmd
        eval $cmd
    else # phase error
        echo "phase"
        cmd="python BEAM_error_difference.py 8192 2048 10.0 1.0 3.0 y phase $element"
        echo $cmd
        eval $cmd
    fi
    echo "Done"
    echo
done

echo




# python BEAM_test_new.py 12288 1024 10 1.0 n 1.0
# python BEAM_test_new.py 12288 1024 10 1.0 n 2.0
# python BEAM_test_new.py 12288 1024 10 1.0 n 3.0
# python BEAM_test_new.py 12288 1024 10 1.0 n 4.0
# python BEAM_test_compare.py 1024 3.0