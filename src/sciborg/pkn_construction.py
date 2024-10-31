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
    from .pyBRAvo.src import  pyBravo
    import shutil

except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()


def output_dir_for_pybravo(dir_:str) -> str:
    return f'{dir_}/{os.path.basename(dir_)}'


def format_args(args:dict) -> dict:
    formated_args = dict()
    formated_args['f'] = args['inputs_genes_file']
    formated_args['reg'] = True if args['reconstruction_type']=='regulation' else False
    formated_args['sig'] = True if args['reconstruction_type']=='signaling' else False
    formated_args['md'] = int(args['max_depth'])
    formated_args['endpoint'] = args['endpoint']
    formated_args['excl'] = str()
    formated_args['excl'] = args['exclude_databases']
    formated_args['s'] = True if args['extend_with_synonyms'] else False
    formated_args['su'] = True if args['extend_with_rna_protein_suffixes'] else False
    formated_args['c'] = True if args['decompose_complexes'] else False
    formated_args['unk'] = True if args['unknown'] else False
    formated_args['fast'] = True if args['fast'] else False
    formated_args['o'] = output_dir_for_pybravo(args['output_dir'])
    # None used pyBRAvo params
    formated_args['sigd'] = False
    formated_args['v'] = False
    formated_args['i'] = None
    formated_args['incl'] = None
    return formated_args





def run_pybravo(args):
    

    # Copy the input genes file
    src_ = args['inputs_genes_file']
    target = f"{args['output_dir']}/input_genes_list.txt"
    shutil.copy(src_, target)

    # Format args for pyBRAvo
    args = format_args(args)
    args = Namespace(**args)
    pyBravo.main(args)
