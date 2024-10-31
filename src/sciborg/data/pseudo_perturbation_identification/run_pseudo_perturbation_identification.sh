#! /bin/sh
make_copy () {
    i=1
    to_copy=${out_dir}/pseudo_perturbation_answer_sets.txt
    while  [[ ! $(tail -n 1 $to_copy) =~ ^CPU.* ]]
    do
        sleep $(($1));
        cp ${to_copy} ${out_dir}/timelaps/$(($i*$1)).out
        ((i++))

    done
}

k=$1
i=$2
out_dir=$3
timeout=$4
expression_instance_filename=$5
input_instance_filename=$6
parents_instance_filename=$7
intermediates_instance_filename=$8
problem_instance=$9
time_point=${10}


mkdir ${out_dir}/timelaps

make_copy $time_point & clingo --const k=${k} --const i=${i} --time-limit=${timeout} ${expression_instance_filename} ${input_instance_filename} ${parents_instance_filename} ${intermediates_instance_filename} ${problem_instance} > ${out_dir}/pseudo_perturbation_answer_sets.txt