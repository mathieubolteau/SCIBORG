#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    import csv
    import pandas as pd
    import numpy as np 
    import re
    import json
    from argparse import Namespace
    import pydotplus
    import os
    import networkx as nx
    import matplotlib.pyplot as plt

    from .utils import load_data, read_file, get_cell_class
    from caspo.console.handlers import learn_handler, classify_handler, visualize_handler


except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()


def find_last_index(search_list, starting_item):
    starting_item = '^'+starting_item
    indices = [ i for i, word in enumerate(search_list) if re.search(starting_item, word) ]
    return indices[-1]

def get_value(atom):
    return re.findall('\"([A-Za-z0-9\.\-]*)\"', atom)


# def get_value(atom):
#     print('atom:', atom)
#     print(re.findall('\"([A-Z0-9\.]*)\"', atom))
#     return re.findall('([0-9]+)', atom)


def load_answer_set(filename: str) -> set:
    sel_genes = list()
    affinities = list()
    file_list = read_file(filename)
    index = find_last_index(file_list, 'Optimization:')-1 # -1 to have the answer-set of the last optimization
    answer_set_raw = file_list[index]
    answer_set_list = answer_set_raw.split(' ')
    for atom in answer_set_list:
        if re.search('^selgene', atom):
            encoded_gene = get_value(atom)[0]
            gene = encoded_gene
            sel_genes.append(gene)
        elif re.search('^affinity', atom):
            cur_cells = get_value(atom)
            encoded_cell_i = cur_cells[0]
            cell_i = encoded_cell_i
            encoded_cell_j = cur_cells[1]
            cell_j = encoded_cell_j
            affinities.append([cell_i, cell_j])
            

    return sel_genes, affinities

def get_classes_cells(matrix, classes):
    # Not necessary, see config dict
    classes_cells = dict()
    for class_ in classes:
        classes_cells[class_] = matrix[matrix['clusterUmap']==class_].index.tolist()
    return classes_cells


def get_expr_vector(cell:str, genes:list, matrix:pd.DataFrame)-> list:
    expr_vector = list()
    for gene in genes:
        if gene in matrix.columns.to_list(): # TO REMOVE
            expr_vector.append(matrix.at[cell, gene])
    return expr_vector

def get_redundancies(studied_cell:str, selgenes:list, class_:str, matrix:pd.DataFrame)-> list:
    redundancies = list()
    redundancies.append(studied_cell)
    cells = matrix.index.tolist()
    studied_cell_vector = get_expr_vector(studied_cell, selgenes, matrix)
    for cell in cells: 
        if cell != studied_cell:

            vector = get_expr_vector(cell, selgenes, matrix)
            if vector == studied_cell_vector:
                redundancies.append(cell)
    
    return redundancies

def redundancies_calculation(affinities:dict, selgenes:list, class_names:list, matrix:pd.DataFrame)-> dict:
    redundancy_cells = {
        class_names[0]: {},
        class_names[1]: {}
    }
    redundancy_vectors = {}
    c1_class_name = class_names[0]
    c1_submatrix = matrix[matrix['clusterUmap']==c1_class_name]

    c2_class_name = class_names[1]
    c2_submatrix = matrix[matrix['clusterUmap']==c2_class_name]
    for i in range(len(affinities)):
        c1 = affinities[i][0]
        c2 = affinities[i][1]
        redundancy_vectors[i] = list()
        c1_redundancies = get_redundancies(c1, selgenes, c1_class_name, c1_submatrix)
        redundancy_cells[c1_class_name][i] = c1_redundancies
        c2_redundancies = get_redundancies(c2, selgenes, c2_class_name, c2_submatrix)
        redundancy_cells[c2_class_name][i] = c2_redundancies
        
        for c1_red in c1_redundancies:
            for c2_red in c2_redundancies:
                redundancy_vectors[i].append((c1_red, c2_red))
    
    sum_ = 0
    for values in redundancy_cells[c1_class_name].values():
        sum_ += len(values)
    redundancy_cells[c1_class_name]['nb'] = sum_
    sum_ = 0
    for values in redundancy_cells[c2_class_name].values():
        sum_ += len(values)
    redundancy_cells[c2_class_name]['nb'] = sum_
    return redundancy_cells, redundancy_vectors



def readouts_diff_calculation(aff_i, aff_j, readouts, matrix) -> float:
    aff_i_vector = get_expr_vector(aff_i,readouts, matrix)
    aff_j_vector = get_expr_vector(aff_j,readouts, matrix)
    diff = 0
    for idx in range(len(aff_i_vector)):
        diff += abs(aff_i_vector[idx]-aff_j_vector[idx])
    return diff

def readouts_maximization(redundancy_vector:dict, readouts:list, matrix:pd.DataFrame) -> list:
    perturbations = list()
    for aff in range(len(redundancy_vector)):
        max_pair = 0
        max_idx = 0
        for idx in range(len(redundancy_vector[aff])):
            aff_i = redundancy_vector[aff][idx][0]
            aff_j = redundancy_vector[aff][idx][1]
            readout_diff = readouts_diff_calculation(aff_i, aff_j, readouts, matrix)
            if readout_diff > max_pair:
                max_pair = readout_diff
                max_idx = idx 
        perturbations.append(redundancy_vector[aff][max_idx])
    return perturbations


def create_midas_file(expr_data: pd.DataFrame, readout_genes: str, input_genes: list, cells: list, out_dir: str, class_name: str, intermediate_genes:list) :

    inputs_dict = dict()
    DA_readouts_dict = dict()
    DV_readouts_dict = dict()

    class_len = len(cells)

    inputs_dict[f'TR:{class_name}:CellLine'] = [1 for x in range(class_len*2)]

    for ii_gene in input_genes:
        current_column = list()
        for cell in cells:
            expr_value = expr_data.at[cell, ii_gene]
            current_column.append(expr_value)
        inputs_dict[f'TR:{ii_gene}'] = current_column+current_column
    for ii_gene in intermediate_genes:
        current_column = list()
        for cell in cells:
            expr_value = expr_data.at[cell, ii_gene]
            expr_value = 1 if expr_value==0 else 0
            current_column.append(expr_value)
        inputs_dict[f'TR:{ii_gene}'] = current_column+current_column

    for r_gene in readout_genes:
        DA_readouts_dict[f'DA:{r_gene}'] = [0 for x in range(class_len)]+[10 for x in range(class_len)]

        current_column = list()
        for cell in cells:
            current_column.append(expr_data.at[cell, r_gene])

        DV_readouts_dict[f'DV:{r_gene}'] = [0 for x in range(class_len)]+current_column

    inputs_df = pd.DataFrame(inputs_dict)
    DA_readouts_df = pd.DataFrame(DA_readouts_dict)
    DV_readouts_df = pd.DataFrame(DV_readouts_dict)

    # concatenate the 3 dataframe
    midas_df = pd.concat([inputs_df, DA_readouts_df, DV_readouts_df], axis=1)
    midas_df.to_csv(f'{out_dir}/{class_name}_midas.csv', index=False)
 

def create_setup_file(ii_genes:list, readouts_genes:list, input_genes:list, intermediates_genes:list, out_dir:str):
    setup = {
        'stimuli' : list(),
        'inhibitors' : list(),
        'readouts' : list()
    }
    setup['readouts'] = readouts_genes
    for gene in ii_genes:
        if gene in input_genes:
            setup['stimuli'].append(gene)
        elif gene in intermediates_genes:
            setup['inhibitors'].append(gene)
    json.dump(setup, open(f"{out_dir}/setup.json", 'w'))

def transform_sif_file(out_dir:str):
    with open(f'{out_dir}/reduced_pkn.sif') as f:
        content = f.read()
    to_replace = {
        'PART_OF' : '1',
        'ACTIVATION' : '1',
        'INHIBITION' : '-1',
        '\"' : '',
        ' ': '',
        '.*\sUNKNOWN\s.*\n':''
    }
    for key,value in to_replace.items():
        # content = content.replace(key, value)
        content = re.sub(key, value, content)
    with open(f'{out_dir}/reduced_transformed_pkn.sif', 'w') as f:
        f.write(content)


def learn_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = f"{args['output_dir']}"
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/reduced_transformed_pkn.sif'
    formatted_args['midas'] = f'{out_dir}/{class_}_midas.csv'
    formatted_args['time'] = 10
    formatted_args['fit'] = args['fit'] if args['fit'] else 0
    formatted_args['size'] = args['size'] if args['size'] else 0
    formatted_args['length'] = args['length'] if args['length'] else 0
    formatted_args['optimum'] = args['optimum'] if args['optimum'] else 0
    # Non-used params - Caspo default value
    formatted_args['discretization'] = 'round'  
    formatted_args['factor'] = 100             
    formatted_args['threads'] = 1              
    formatted_args['conf'] = 'many'            
    return Namespace(**formatted_args)


def classify_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = f"{args['output_dir']}"
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/reduced_transformed_pkn.sif'
    formatted_args['setup'] = f'{out_dir}/setup.json'
    formatted_args['networks'] = f'{out_dir}/{class_}/networks.csv'
    formatted_args['midas'] = [f'{out_dir}/{class_}_midas.csv', 10]
    # Non-used params - Caspo default value   
    formatted_args['threads'] = 1              
    formatted_args['conf'] = 'many'         
    return Namespace(**formatted_args)


def visualize_format_args(args:dict, class_:str)-> Namespace:
    formatted_args = dict()
    out_dir = f"{args['output_dir']}"
    formatted_args['out'] = f'{out_dir}/{class_}'
    formatted_args['pkn'] = f'{out_dir}/reduced_transformed_pkn.sif'
    formatted_args['setup'] = f'{out_dir}/setup.json'
    formatted_args['midas'] = [f'{out_dir}/{class_}_midas.csv', 10]
    formatted_args['networks'] = f'{out_dir}/{class_}/networks.csv'
    formatted_args['sample'] = args['sample'] if args['sample'] else -1
    # Non-used params  
    formatted_args['stats_networks'] = None              
    formatted_args['stats_strategies'] = None           
    formatted_args['strategies'] = None 
    formatted_args['designs'] = None 
    formatted_args['behaviors'] = None 
    formatted_args['predictions'] = None 
    return Namespace(**formatted_args)



def dot_to_pdf(config:dict, class_:str):
    out_dir = f'{config["output_dir"]}/{class_}'
    # TODO Handle sample = [0,N]
    # if config['sample'] != '-1' :
    #     to_convert = ['networks-union', 'pkn', 'pkn-zip']
    # else :
    #     to_convert = ['pkn', 'pkn-zip']
    to_convert = ['networks-union', 'pkn', 'pkn-zip']
    for file in to_convert:
        in_path = f'{out_dir}/{file}.dot'
        out_path = f'{out_dir}/{file}.pdf'
        graph = pydotplus.graphviz.graph_from_dot_file(in_path)
        graph.write_pdf(out_path)

def get_timelaps_redundancies(filename:str, classes:list, expr_data:pd.DataFrame):
    sel_genes, affinities = load_answer_set(filename)

    redundancies_cells, _ = redundancies_calculation(affinities, sel_genes, classes, expr_data)
    return redundancies_cells[classes[0]]['nb'], redundancies_cells[classes[1]]['nb']


def timelaps_analyze(directory:str, classes:list, expr_data:pd.DataFrame, out_dir:str):
    steps = os.listdir(directory)
    steps = [x  for x in steps if x.endswith(".out")]
    steps = sorted([int(step_file.split(".")[0]) for step_file in steps])

    X = list()
    Y_affinities = list()
    Y_redundancies_class1 = list()
    Y_redundancies_class2 = list()
    
    for step in steps:
        step_file = f"{directory}/{step}.out"
        with open(step_file) as f:
            lines = f.readlines()
        opti = False
        last_optimization_line = None
        for line in reversed(lines):
            if line.startswith("Optimization:"):
                last_optimization_line = line
                break
        if last_optimization_line:
            opti = int(last_optimization_line.split("Optimization:")[-1].strip())
            index = lines.index(last_optimization_line)
            answer_set = lines[index - 1].strip()
            X.append(step)
            Y_affinities.append(opti*-1)
            red_class1, red_class2 = get_timelaps_redundancies(step_file, classes, expr_data)
            Y_redundancies_class1.append(red_class1)
            Y_redundancies_class2.append(red_class2)
    

    # Create Plot
    fig, ax1 = plt.subplots() 
    
    ax1.set_xlabel('Time') 
    ax1.set_ylabel('Number of affinity') 
    line1, = ax1.plot(X, Y_affinities, color = '#66c2a5', marker='o', label='Number of affinity') 
    ax1.tick_params(axis ='y') 

    
    # Adding Twin Axes

    ax2 = ax1.twinx() 
    ax1.set_xticks(X)
    ax2.set_ylabel('Number of redundancy') 
    line2, = ax2.plot(X, Y_redundancies_class1, color = '#8da0cb', marker='o', label='CLASS 1') 
    line3, = ax2.plot(X, Y_redundancies_class2, color = '#8da0ee', marker='o', label='CLASS 2') 
    ax2.tick_params(axis ='y') 
    ax1.legend(handles=[line1, line2, line3])
    plt.savefig(f'{out_dir}/redundancies_two_classes.png')

    # Create Plot
    Y_redundancies = [a + b for a, b in zip(Y_redundancies_class1, Y_redundancies_class2)]
    fig, ax1 = plt.subplots() 
    
    ax1.set_xlabel('Time') 
    ax1.set_ylabel('Number of affinity') 
    line1, = ax1.plot(X, Y_affinities, color = '#66c2a5', marker='o', label='Number of affinity') 
    ax1.tick_params(axis ='y')
    
    # Adding Twin Axes

    ax2 = ax1.twinx() 
    ax1.set_xticks(X)
    ax2.set_ylabel('Number of redundancy') 
    line2, = ax2.plot(X, Y_redundancies, color = '#8da0cb', marker='o', label='redundancy') 
    ax2.tick_params(axis ='y') 
    ax1.legend(handles=[line1, line2])
    plt.savefig(f'{out_dir}/redundancies_merged_classes.png')





def run_bn_inference(config):
    # out_dir = config['output_dir']
    classes = config['class_types']
    # matrix_filename = f"{out_dir}/bin_reduced_matrix.csv"

    # expr_data = load_data(matrix_filename, index_name='Name')
    # answer_set_filename = f"{out_dir}/pseudo_perturbation_answer_sets.txt"
    # sel_genes, affinities = load_answer_set(answer_set_filename)

    # input_genes_in_the_matrix = read_file(f"{out_dir}/no_predecessors_in_the_matrix.txt")
    # intermediate_genes_in_the_matrix = read_file(f"{out_dir}/intermediates_in_the_matrix.txt")
    # readout_genes = read_file(f"{out_dir}/no_successors_in_the_matrix.txt")

    # # Exclude gene not in the matrix
    # input_genes = list(set(sel_genes).intersection(input_genes_in_the_matrix))
    # intermediate_genes = list(set(sel_genes).intersection(intermediate_genes_in_the_matrix))
    # input_genes.sort()
    # intermediate_genes.sort()
    # readout_genes.sort()
    
    # # Redundancies calculation 
    # classes = list()
    # for c in affinities[0]:
    #     classes.append(get_cell_class(c, expr_data))
    # classes_cells = get_classes_cells(matrix=expr_data, classes=classes)

    # # Maximization of readouts difference
    # json.dump(affinities, open(f"{out_dir}/affinities.json",'w'))
    
    # redondancy_cells, redondancy_vectors = redundancies_calculation(affinities, sel_genes, classes, expr_data)
    # json.dump(redondancy_cells, open(f"{out_dir}/redondancy_cells.json",'w'))
    # json.dump(redondancy_vectors, open(f"{out_dir}/redondancy_vectors.json",'w'))
    # perturbations = readouts_maximization(redondancy_vectors, readout_genes, expr_data)
    # json.dump(perturbations, open(f"{out_dir}/perturbations.json",'w'))

    # # Preprocessing step
    # for i  in range(len(classes)):
    #     class_ = classes[i]
    #     if not os.path.exists(f'{out_dir}/{class_}'): os.makedirs(f'{out_dir}/{class_}')
    #     curent_cells = [item[i] for item in perturbations] 
    #     create_midas_file(input_genes=input_genes, expr_data=expr_data, cells=curent_cells,
    #                     readout_genes=readout_genes, out_dir=out_dir, class_name=class_, intermediate_genes=intermediate_genes)
    
    # transform_sif_file(out_dir)
    # create_setup_file(sel_genes, readout_genes, input_genes, intermediate_genes, out_dir)

    for class_ in classes:
        print(f"--> Class: {class_}")
        # Learning step
        learn_args = learn_format_args(config, class_)
        learn_handler(learn_args)

        # Classify step
        # classify_args = classify_format_args(config, class_)
        # classify_handler(classify_args)

        # Visualize step
        visualize_args = visualize_format_args(config, class_)
        visualize_handler(visualize_args)
        dot_to_pdf(config, class_)

    