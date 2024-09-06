#!/bin/bash

# This bash file must be run from the directory where it is stored

run_simulation() {
    local reps=$1
    local d_flag=$2
    local t=$3
    local r_block=$4
    local cutloc=$5
    local inst=$6
    local q=$7
    local circtype=$8
    local circname_nocut=$9
    local circname_cut=${10}
    local num_amps=${11}
    local k_values=${12}  # Added for -k flag values

    for (( i=0; i<reps; i++ )); do
        # Run example with cut but no block
        local filename0="./logs_boixo/times_noblock_${circtype}_inst_${inst}_cutloc${cutloc}_rep${i}.tlog"
        local amps0="./logs_boixo/amps_noblock_${circtype}_inst_${inst}_cutloc${cutloc}_rep${i}.log"
        echo "Starting Cut (no block)..."
        timeout 3600 ../../apps/qsimh_amplitudes.x -c "${circname_nocut}" -d "${d_flag}" -k "${k_values}" -t "${t}" -w 0 -p 0 -r "${r_block}" -i "../../circuits/bitstrings_q${q}" -o "$amps0" -v 1 > "$filename0"
        if [ $? -eq 124 ]; then
            echo "Cut (no block) command timed out"
            break
        else
            echo "Completed: $filename0"
        fi
    done

    for (( i=0; i<reps; i++ )); do
        # Run example with cut (block)
        local filename1="./logs_boixo/times_block_${circtype}_inst_${inst}_cutloc${cutloc}_rep${i}.tlog"
        local amps1="./logs_boixo/amps_block_${circtype}_inst_${inst}_cutloc${cutloc}_rep${i}.log"
        echo "Starting Cut (block)..."
        ../../apps/qsimh_amplitudes.x -c "${circname_cut}" -d "${d_flag}" -k "${k_values}" -t "${t}" -w 0 -p 0 -r "${r_block}" -i "../../circuits/bitstrings_q${q}" -o "$amps1" -v 1 > "$filename1"
        echo "Completed: $filename1"
    done

    for (( i=0; i<reps; i++ )); do
        # Run example without cut
        local filename2="./logs_boixo/times_nocut_${circtype}_inst_${inst}_rep${i}.tlog"
        local amps2="./logs_boixo/amps_nocut_${circtype}_inst_${inst}_rep${i}.log"
        echo "Starting no cut..."
        ../../apps/qsim_amplitudes.x -c "${circname_nocut}" -d "${d_flag}" -i "../../circuits/bitstrings_q${q}" -o "$amps2" -v 2 -t "${t}" > "$filename2"
        echo "Completed: $filename2"
    done
}

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=11
inst="4x5_10_0"  # "5x5_22_9"
q=20             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=11
inst="4x5_11_0"  # "5x5_22_9"
q=20             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=11
inst="4x5_12_0"  # "5x5_22_9"
q=20             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=11
inst="5x5_10_0"  # "5x5_22_9"
q=25             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=11
inst="5x5_11_0"  # "5x5_22_9"
q=25             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=11
inst="5x5_12_0"  # "5x5_22_9"
q=25             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=13
inst="5x6_10_0"  # "5x5_22_9"
q=30             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=13
inst="5x6_11_0"  # "5x5_22_9"
q=30             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

# Parameters
reps=5
d_flag=2000
t=16
r_block=5
cutloc=13
inst="5x6_12_0"  # "5x5_22_9"
q=30             # must fit inst
circtype="is_v1" # "is_v1" or "cz_v2"
num_amps=1000000 # do not forget to alter num_amps

# -k values
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13"  # Customizable input for the -k flag

# Circuit names
circname_nocut="./${circtype}/inst_${inst}.txt"
circname_cut="./${circtype}/inst_${inst}_block_cutloc${cutloc}.txt"

# Run the simulation with the specified parameters
run_simulation $reps $d_flag $t $r_block $cutloc $inst $q $circtype $circname_nocut $circname_cut $num_amps "$k_values"

echo "Finished: $inst"

#------------------------------------------------

