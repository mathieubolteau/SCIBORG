#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    import csv
    import pandas as pd
    import matplotlib.pyplot as plt
    from scipy import stats
    import json
    import os
    from sys import argv
    from argparse import Namespace
    import configparser
    # from .pyBRAvo.src import pyBravo
    import shutil
    from math import ceil
    import pkg_resources


    from . import preprocessing


except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()




def run_pseudo_perturbation_identification(config):
    out_dir = config['output_dir']
    k = config['k']
    i = ceil(k*config['gene_inputs_pourcentage']/100)
    expression_instance_filename = f'{out_dir}/expression_instance.lp'
    input_instance_filename = f'{out_dir}/inputs_instance.lp'
    intermediates_instance_filename = f'{out_dir}/intermediates_instance.lp'
    ancestors_instance_filename = f'{out_dir}/ancestors_instance.lp'
    problem_instance = pkg_resources.resource_filename(__name__, 'data/pseudo_perturbation_identification/problem.lp')
    timeout = config['timeout'] if config['timeout'] != '' else 0
    time_point = config['time_point']
    script = pkg_resources.resource_filename(__name__, 'data/pseudo_perturbation_identification/run_pseudo_perturbation_identification.sh')

    # Preprocessing step, considering that data from PKN construction are available
    # TODO: check that the PKN construction step is done before this step.
    preprocessing.run_preprocessing(config)


    # Run the pseudo-perturbation identification program 
    # cmd = f'clingo --const k={k} --const i={i} --time-limit={timeout} {expression_instance_filename} {input_instance_filename} {ancestors_instance_filename} {intermediates_instance_filename} {problem_instance} > {out_dir}/pseudo_perturbation_answer_sets.txt'
    cmd = f'bash {script} {k} {i} {out_dir} {timeout} {expression_instance_filename} {input_instance_filename} {ancestors_instance_filename} {intermediates_instance_filename} {problem_instance} {time_point}'
    print(f'CMD:  {cmd}')
    clingo_output = os.popen(cmd).read()
    print(clingo_output)