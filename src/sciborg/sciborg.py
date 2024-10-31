#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    # Standard library imports
    import argparse
    import configparser
    import logging
    import pandas as pd
    from sys import argv
    import json
    import os

    # Local imports
    from . import utils
    from . import pkn_construction
    # from . import preprocessing
    from . import pkn_analyze
    from . import asp_analyze
    from . import pseudo_perturbation_identification
    from . import pseudo_observation_diff_maxi
    from . import bn_inference

except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()

def config_file_parser(configuration_file:str) -> configparser.ConfigParser:
    config = configparser.ConfigParser(allow_no_value=True)
    config.read_file(open(configuration_file))
    return config


def check_config_file(config_parser: configparser.ConfigParser) -> dict():
    config = dict()
    try : 
        # If empty, default value = 'out'
        config['output_dir'] = config_parser.get('PATHS', 'output_directory_path') if config_parser.get('PATHS', 'output_directory_path') != '' else 'out'
        
        
        config['matrix_path'] = config_parser.get('DATA','matrix_path')
        if config['matrix_path'] == '':
            raise ValueError('Configuration file -- [DATA] : Expression matrix path if mandatory. See configuration file.')
        config['annotation_len'] = config_parser.getint('DATA','annotation_len')-1 # -1 to deal with the index start of dataframe equal to 0
        if config['annotation_len'] < 0:
            raise ValueError('Configuration file -- [DATA] : Incorrect annotation len. See configuration file.')
        config['class_types'] = config_parser.get('DATA', 'class_types').split(',') 
        if len(config['class_types']) < 2 or config['class_types']=='': 
            raise ValueError('Configuration file -- [DATA] : 2 classes are required. See configuration file.')

        config['inputs_genes_file'] = config_parser.get('PKN_CONSTRUCTION', 'inputs_genes_file')
        config['reconstruction_type'] = config_parser.get('PKN_CONSTRUCTION', 'reconstruction_type').lower()
        if config['reconstruction_type'] not in ['regulation', 'signalisation']:
            raise ValueError('Configuration file -- [PKN_CONSTRUCTION] : Incorrect reconstruction type. See configuration file.')
        config['max_depth'] = config_parser.getint('PKN_CONSTRUCTION', 'max_depth')
        if config['max_depth'] < 0:
            raise ValueError('Configuration file -- [PKN_CONSTRUCTION] : Incorrect max depth. See configuration file.')
        
        config['extend_with_synonyms'] = config_parser.getboolean('PKN_CONSTRUCTION', 'extend_with_synonyms')
        config['extend_with_rna_protein_suffixes'] = config_parser.getboolean('PKN_CONSTRUCTION', 'extend_with_rna_protein_suffixes')
        config['decompose_complexes'] = config_parser.getboolean('PKN_CONSTRUCTION', 'decompose_complexes')
        config['fast'] = config_parser.getboolean('PKN_CONSTRUCTION', 'fast')
        config['unknown'] = config_parser.getboolean('PKN_CONSTRUCTION', 'unknown')
        config['fast'] = config_parser.getboolean('PKN_CONSTRUCTION', 'fast')
        config['endpoint'] = config_parser.get('PKN_CONSTRUCTION', 'endpoint')
        config['exclude_databases'] = config_parser.get('PKN_CONSTRUCTION', 'exclude_databases').split(',')
        config['exclude_databases'] = [x.lower() for x in config['exclude_databases']]


        config['k'] = config_parser.getint('PSEUDO_PERTURBATIONS', 'k')
        if config['k'] <= 0:
            raise ValueError('Configuration file -- [PSEUDO_PERTURBATIONS] : Incorrect k value. See configuration file.')
        config['gene_inputs_pourcentage'] = config_parser.getint('PSEUDO_PERTURBATIONS', 'gene_inputs_pourcentage')
        if config['gene_inputs_pourcentage'] < 0:
            raise ValueError('Configuration file -- [PSEUDO_PERTURBATIONS] : Incorrect gene inputs pourcentages. See configuration file.')
        config['timeout'] = config_parser.getint('PSEUDO_PERTURBATIONS', 'timeout')
        if config['timeout'] < 0:
            raise ValueError('Configuration file -- [PSEUDO_PERTURBATIONS] : Incorrect timeout. See configuration file.')
        config['time_point'] = config_parser.getint('PSEUDO_PERTURBATIONS', 'time_point')
        if config['time_point'] < 0:
            raise ValueError('Configuration file -- [PSEUDO_PERTURBATIONS] : Incorrect time_point. See configuration file.')

        config['fit'] = config_parser.getfloat('BN_INFERENCE', 'fit')
        if config['fit'] < 0:
            raise ValueError('Configuration file -- [BN_INFERENCE] : Incorrect fit value. See configuration file.')
        config['size'] = config_parser.getint('BN_INFERENCE', 'size')
        if config['size'] < 0:
            raise ValueError('Configuration file -- [BN_INFERENCE] : Incorrect size value. See configuration file.')
        config['length'] = config_parser.getint('BN_INFERENCE', 'length')
        if config['length'] < 0:
            raise ValueError('Configuration file -- [BN_INFERENCE] : Incorrect length value. See configuration file.')
        config['optimum'] = config_parser.getint('BN_INFERENCE', 'optimum')
        if config['optimum'] < 0:
            raise ValueError('Configuration file -- [BN_INFERENCE] : Incorrect optimum value. See configuration file.')
        config['sample'] = config_parser.getint('BN_INFERENCE', 'sample')
        if config['sample'] < -1:
            raise ValueError('Configuration file -- [BN_INFERENCE] : Incorrect sample value. See configuration file.')
        return config

    except configparser.NoOptionError as E:
        print(E)
        print("An option is missing in the configuration file")
        print("Please report to the descriptions in the configuration file\n")
        exit()
    except configparser.NoSectionError as E:
        print(E)
        print("An section is missing in the configuration file")
        print("Please report to the descriptions in the configuration file\n")
        exit()
    except ValueError as E:
        print(E)
        print("One of the value in the configuration file is not in the correct format")
        print("Please report to the descriptions in the configuration file\n")
        exit()


def run():
    sciborg_parser = argparse.ArgumentParser()
    sciborg_parser.add_argument('configuration_file', 
                                 help='configuration file path')
    sciborg_parser.add_argument('-pkn', '--pkn-construction', action='store_true', required=False,
                                help='run the PKN construction')
    # sciborg_parser.add_argument('--preprocessing', action='store_true', required=False,
    #                             help='run the preprocessing step')
    sciborg_parser.add_argument('-pp', '--pseudo-perturbation-identification', action='store_true', required=False,
                                help='run the pseudo-perturbation identification')
    sciborg_parser.add_argument('-po', '--pseudo-observation-diff-maxi', action='store_true', required=False,
                                help='run the pseudo-observation difference maximization')
    sciborg_parser.add_argument('-bn', '--bn-inference', action='store_true', required=False,
                                help='run the BN inference')
    # sciborg_parser.add_argument('-c', '--cell-classifier', action='store_true', required=False,
    #                             help='NOT AVAILABLE -- run the cell classifier [Work in progress]')
    

    args = sciborg_parser.parse_args()
    config_parser = config_file_parser(args.configuration_file)
    config = check_config_file(config_parser)
    out_dir_path = config['output_dir']
    out_dir_name = os.path.basename(os.path.normpath(out_dir_path))
    
    # Create output directory
    if not os.path.isdir(out_dir_path):
        os.mkdir(out_dir_path)

    # PKN construction
    if args.pkn_construction:
        pkn_construction.run_pybravo(args=config)
        sif_file = f"{out_dir_path}/{out_dir_name}-unified.sif"
        modified_sif_file = f"{out_dir_path}/reduced_pkn.sif"
        gene_expr_mtx_file = config['matrix_path']
        input_genes_file = config['inputs_genes_file']
        annotation_len = config['annotation_len']
        pkn_analyze.run_pkn_analyze(sif_file, out_dir_path, gene_expr_mtx_file, input_genes_file, annotation_len)
        asp_analyze.run_asp_analyze(modified_sif_file, gene_expr_mtx_file, out_dir_path)

    # Preprocessing step, considering that data from PKN construction are available
    # TODO: check that the PKN construction step is done before this step.
    # if args.preprocessing:
    #     preprocessing.run_preprocessing(config)
    
    # Pseudo-perturbation identification step
    if args.pseudo_perturbation_identification:
        pseudo_perturbation_identification.run_pseudo_perturbation_identification(config)

    # Pseudo-perturbation identification step
    if args.pseudo_observation_diff_maxi:
        pseudo_observation_diff_maxi.run_pseudo_observation_diff_maxi(config)


    # BNs inference step
    if args.bn_inference:
        bn_inference.run_bn_inference(config)

        




if __name__ == '__main__':
    run()