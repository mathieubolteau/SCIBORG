import pandas as pd
from sys import argv
import csv
import clyngor
import json
import os
import pkg_resources

def sif2lp(path):
    encoding = ''
    df = pd.read_csv(path, delimiter='\t', names=['source', 'sign', 'target']).drop_duplicates()
    cpt=1
    for _, source, sign, target in df.itertuples():
        encoding += f'reaction("r{cpt}","{sign}"). reactant("{source}", "r{cpt}"). product("{target}", "r{cpt}").\n'
        cpt += 1

    return encoding

def to_file(str, outname):
    with open(outname, 'w') as f:
        f.write(str)


def get_no_predecessors(encoding):
    instance_path = pkg_resources.resource_filename(__name__, 'data/pkn_construction/get_no_predecessors.lp')
    answers = clyngor.solve(instance_path, inline=encoding)
    for answer in answers.by_predicate.first_arg_only:
        no_predecessors = list(answer['no_predecessor'])
        no_predecessors = [x.replace('"','') for x in no_predecessors]
        return no_predecessors

def get_no_successors(encoding):
    instance_path = pkg_resources.resource_filename(__name__, 'data/pkn_construction/get_no_successors.lp')
    answers = clyngor.solve(instance_path, inline=encoding)
    for answer in answers.by_predicate.first_arg_only:
        no_successors = list(answer['no_successor'])
        no_successors = [x.replace('"','') for x in no_successors]
        return no_successors

def get_complexes(encoding):


    # TODO : fix the import of 'get_complexes.lp' file.
    instance_path = pkg_resources.resource_filename(__name__, 'data/pkn_construction/get_complexes.lp')
    
    answers = clyngor.solve(instance_path, inline=encoding)
    for answer in answers.by_predicate.first_arg_only:
        complexes = list(answer['complexes'])
        complexes = [x.replace('"','') for x in complexes]
        return complexes

def get_nodes(encoding):
    instance_path = pkg_resources.resource_filename(__name__, 'data/pkn_construction/get_nodes.lp')
    answers = clyngor.solve(instance_path, inline=encoding)
    for answer in answers.by_predicate.first_arg_only:
        nodes = list(answer['node'])
        nodes = [x.replace('"','') for x in nodes]
        return nodes

def count_nb_edges(encoding):
    # Each line of the encoding program correspond to a reaction, so an edge
    return len(encoding.split('\n'))

def count_edges_type(encoding):
    instance_path = pkg_resources.resource_filename(__name__, 'data/pkn_construction/count_edges_type.lp')
    answers = clyngor.solve(instance_path, inline=encoding)
    for answer in answers.by_predicate:
        edges_type = dict()
        for tuple_ in answer['nb']:
            edges_type[tuple_[0].replace('"','')] = tuple_[1]
        return edges_type

def count_nb_lines(file):
    return len(open(file, 'r').readlines())

def compute_stats(dir, encoding, nodes, no_predecessors, no_successors):
    stats = dict()
    stats['nb_nodes'] = len(nodes)
    stats['nb_edges'] = count_nb_edges(encoding)
    stats['edges'] = count_edges_type(encoding)
    stats['nb_input_genes'] = count_nb_lines(dir+'/input_genes_list.txt')
    stats['nb_input_genes_not_in_the_matrix'] = count_nb_lines(dir+'/analyze/input_genes_not_in_the_pkn.txt')
    stats['nb_genes_not_in_expr_mtx'] = count_nb_lines(dir+'/analyze/genes_not_in_expr_mtx.txt')
    stats['nb_genes_almost_in_expr_mtx'] = count_nb_lines(dir+'/analyze/genes_almost_in_expr_mtx.txt')
    stats['nb_nodes_without_predecessors'] = len(no_predecessors)
    stats['nb_noddes_without_successors'] = len(no_successors)
    return stats



def run_asp_analyze(sif_file:str, matrix_filename:str, out_dir:str):
    if not os.path.exists(f'{out_dir}/visu_cys'): os.makedirs(f'{out_dir}/visu_cys')

    encoding = sif2lp(f'{out_dir}/reduced_pkn.sif')
    # perfect_encoding = sif2lp(f'{out_dir}/perfect_matchs/reduced_pkn.sif')



    to_file(encoding, f'{out_dir}/pkn.lp')
    nodes = get_nodes(encoding)
    no_predecessors = get_no_predecessors(encoding)
    no_successors = get_no_successors(encoding)
    complexes = get_complexes(encoding)
    intermediates = list( set(nodes) - set(no_predecessors+no_successors+complexes) )

    matrix_filename = f'{out_dir}/pkn_gene_reduced_expr_mtx.csv'

    matrix_genes = pd.read_csv(matrix_filename).columns

    no_predecessors_in_the_matrix = list( set(no_predecessors).intersection(set(matrix_genes)) )
    no_successors_in_the_matrix = list( set(no_successors).intersection(set(matrix_genes)) )
    intermediates_in_the_matrix = list( set(intermediates).intersection(set(matrix_genes)) )

    stats = compute_stats(out_dir, encoding, nodes, no_predecessors, no_successors)

    with open(f"{out_dir}/nodes.txt", "w") as outfile:
            for item in nodes:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/no_predecessors.txt", "w") as outfile:
            for item in no_predecessors:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/no_successors.txt", "w") as outfile:
            for item in no_successors:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/complexes.txt", "w") as outfile:
            for item in complexes:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/intermediates.txt", "w") as outfile:
            for item in intermediates:
                outfile.write("{}\n".format(item))

    with open(f"{out_dir}/no_predecessors_in_the_matrix.txt", "w") as outfile:
            for item in no_predecessors_in_the_matrix:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/no_successors_in_the_matrix.txt", "w") as outfile:
            for item in no_successors_in_the_matrix:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/intermediates_in_the_matrix.txt", "w") as outfile:
            for item in intermediates_in_the_matrix:
                outfile.write("{}\n".format(item))
    with open(f"{out_dir}/pkn_stats.json", 'w') as f:
        json.dump(stats, f, indent=4, sort_keys=True)


    with open(f"{out_dir}/visu_cys/node_table.csv", 'w') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        header = ['shared name', 'name', 'feature']
        csvwriter.writerow(header)
        for node in nodes:
            if node in no_predecessors: 
                row = [node, node, 'no_predecessors']
            elif node in no_successors : 
                row = [node, node, 'no_successors']
            elif node in complexes:
                row = [node, node,'complexes']
            else: 
                row = [node, node, 'intermediates']
            csvwriter.writerow(row)



if __name__ == '__main__':

    sif_file = argv[1]
    matrix_filename = argv[2]
    out_dir = argv[3]

    encoding = sif2lp(sif_file)


    to_file(encoding, out_dir+'/pkn.lp')
    nodes = get_nodes(encoding)
    no_predecessors = get_no_predecessors(encoding)
    no_successors = get_no_successors(encoding)
    intermediates = list( set(nodes) - set(no_predecessors+no_successors) )

    matrix_genes = pd.read_csv(matrix_filename).columns

    no_predecessors_in_the_matrix = list( set(no_predecessors).intersection(set(matrix_genes)) )
    no_successors_in_the_matrix = list( set(no_successors).intersection(set(matrix_genes)) )
    intermediates_in_the_matrix = list( set(intermediates).intersection(set(matrix_genes)) )

    stats = compute_stats(out_dir, encoding, nodes, no_predecessors, no_successors)

    with open(out_dir+"/nodes.txt", "w") as outfile:
            for item in nodes:
                outfile.write("{}\n".format(item))
    with open(out_dir+"/no_predecessors.txt", "w") as outfile:
            for item in no_predecessors:
                outfile.write("{}\n".format(item))
    with open(out_dir+"/no_successors.txt", "w") as outfile:
            for item in no_successors:
                outfile.write("{}\n".format(item))
    with open(out_dir+"/intermediates.txt", "w") as outfile:
            for item in intermediates:
                outfile.write("{}\n".format(item))

    with open(out_dir+"/no_predecessors_in_the_matrix.txt", "w") as outfile:
            for item in no_predecessors_in_the_matrix:
                outfile.write("{}\n".format(item))
    with open(out_dir+"/no_successors_in_the_matrix.txt", "w") as outfile:
            for item in no_successors_in_the_matrix:
                outfile.write("{}\n".format(item))
    with open(out_dir+"/intermediates_in_the_matrix.txt", "w") as outfile:
            for item in intermediates_in_the_matrix:
                outfile.write("{}\n".format(item))
    with open(out_dir+'/pkn_stats.json', 'w') as f:
        json.dump(stats, f, indent=4, sort_keys=True)


    with open(out_dir+"/visu_cys/node_table.csv", 'w') as f:
        csvwriter = csv.writer(f, delimiter='\t')
        header = ['shared name', 'name', 'feature']
        csvwriter.writerow(header)
        for node in nodes:
            if node in no_predecessors: 
                row = [node, node, 'no_predecessors']
            elif node in no_successors : 
                row = [node, node, 'no_successors']
            else: 
                row = [node, node, 'NA']
            csvwriter.writerow(row)

