#!/bin/bash

run_simulation() {
    local d_flag=$1
    local q=$2
    local t=$3
    local r_block=$4
    local cutloc=$5
    local seed=$6
    local sizes=$7
    local angles=$8
    local pinter=$9
    local pintra=${10}
    local k_values=${11}
    local reps=${12:-2}

    local circname_nocut="grouping_noblock_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes_${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights"
    local circname_cut="grouping_block_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes_${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights"

    for (( i=0; i<reps; i++ )); do
        # Run example without cut
        local filename2="./logs/times_nocut_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights_rep${i}.tlog"
        local amps2="./logs/amps_nocut_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights_rep${i}.log"
        echo "Starting no cut..."
        timeout 3600 ../apps/qsim_amplitudes.x -c ${circname_nocut} -d ${d_flag} -i ../circuits/bitstrings_q${q} -o "$amps2" -v 2 -t ${t} > "$filename2"
        if [ $? -eq 124 ]; then
            echo "No cut command timed out"
            break
        else
            echo "Completed: $filename2"
        fi
    done

    for (( i=0; i<reps; i++ )); do
        # Run example with cut (block)
        local filename1="./logs/times_block_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights_rep${i}.tlog"
        local amps1="./logs/amps_block_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights_rep${i}.log"
        echo "Starting Cut (block)..."
        timeout 3600 ../apps/qsimh_amplitudes.x -c ${circname_cut} -d ${d_flag} -k ${k_values} -t ${t} -w 0 -p 0 -r ${r_block} -i ../circuits/bitstrings_q${q} -o "$amps1" -v 1 > "$filename1"
        if [ $? -eq 124 ]; then
            echo "Cut (block) command timed out"
            break
        else
            echo "Completed: $filename1"
        fi
    done

    for (( i=0; i<reps; i++ )); do
        # Run example with cut but no block
        local filename0="./logs/times_noblock_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights_rep${i}.tlog"
        local amps0="./logs/amps_noblock_q${q}_cutloc${cutloc}_blockgraph_seed${seed}_sizes${sizes}_angles_${angles}_pinter${pinter}_pintra${pintra}_weights_rep${i}.log"
        echo "Starting Cut (no block)..."
        timeout 3600 ../apps/qsimh_amplitudes.x -c ${circname_nocut} -d ${d_flag} -k ${k_values} -t ${t} -w 0 -p 0 -r ${r_block} -i ../circuits/bitstrings_q${q} -o "$amps0" -v 1 > "$filename0"
        if [ $? -eq 124 ]; then
            echo "Cut (no block) command timed out"
            break
        else
            echo "Completed: $filename0"
        fi
    done
}


#---------------------------------------------
d_flag=2000
q=30
t=16
r_block=5
cutloc=14
seed=4
sizes="15_15"
angles="155739_232636"
pinter="010"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps


#---------------------------------------------
d_flag=2000
q=30
t=16
r_block=5
cutloc=14
seed=4
sizes="15_15"
angles="155739_232636"
pinter="015"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps

#---------------------------------------------
d_flag=2000
q=30
t=16
r_block=5
cutloc=14
seed=4
sizes="15_15"
angles="155739_232636"
pinter="017"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps

#---------------------------------------------
d_flag=2000
q=31
t=16
r_block=5
cutloc=14
seed=4
sizes="15_16"
angles="155739_232636"
pinter="010"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps

#---------------------------------------------
d_flag=2000
q=31
t=16
r_block=5
cutloc=14
seed=4
sizes="15_16"
angles="155739_232636"
pinter="015"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps


#---------------------------------------------
d_flag=2000
q=31
t=16
r_block=5
cutloc=14
seed=4
sizes="15_16"
angles="155739_232636"
pinter="017"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps

#---------------------------------------------
d_flag=2000
q=32
t=16
r_block=5
cutloc=15
seed=4
sizes="16_16"
angles="155739_232636"
pinter="010"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps


#---------------------------------------------
d_flag=2000
q=32
t=16
r_block=5
cutloc=15
seed=4
sizes="16_16"
angles="155739_232636"
pinter="011"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps


#---------------------------------------------
d_flag=2000
q=32
t=16
r_block=5
cutloc=15
seed=4
sizes="16_16"
angles="155739_232636"
pinter="012"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps


#---------------------------------------------
d_flag=2000
q=33
t=16
r_block=5
cutloc=15
seed=4
sizes="16_17"
angles="155739_232636"
pinter="009"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps

#---------------------------------------------
d_flag=2000
q=33
t=16
r_block=5
cutloc=15
seed=4
sizes="16_17"
angles="155739_232636"
pinter="010"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps




#---------------------------------------------
d_flag=2000
q=33
t=16
r_block=5
cutloc=15
seed=4
sizes="16_17"
angles="155739_232636"
pinter="011"
pintra="080"
k_values="0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15"
reps=5

# Call the function with the stored parameters
run_simulation $d_flag $q $t $r_block $cutloc $seed $sizes $angles $pinter $pintra "$k_values" $reps
