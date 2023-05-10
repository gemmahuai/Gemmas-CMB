#!/bin/bash

error="${1}" # phase or amplitude

array=(1.0 5.0 10.0 30.0 50.0 100.0)
for element in "${array[@]}"
do
    echo "amplitude = $element"
    if [ "$error" == "amp" ]
    then
        cmd="amp $element"
        echo $cmd
        #eval $cmd
    else
        cmd="phase $element"
        echo $cmd
    fi
done

#echo "1/${#array[@]}"
echo