import csv
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import json
import os
from sys import argv
import networkx as nx 
import re

from .utils import index_syn, init_gene_synonyms_cache, read_file, save_to_file, load_data

def SIF_to_digraph(filename:str) -> nx.DiGraph:
    G = nx.Digraph()
    with open(filename, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for row in csvreader:
            source = row[0]
            edge_type = row[1]
            target = row[2]
            G.add_node(source)
            G.add_node(target)
            G.add_edge(source, target, label=edge_type)
            # raw_data.append(row)
            # if row[0] not in genes_list: genes_list.append(row[0])
            # if row[2] not in genes_list: genes_list.append(row[2])
    return G

def get_genes(sif_file):

    raw_data = list()
    genes_list = list()

    with open(sif_file, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t')
        for row in csvreader:
            raw_data.append(row)
            if row[0] not in genes_list: genes_list.append(row[0])
            if row[2] not in genes_list: genes_list.append(row[2])
    return (raw_data, genes_list)
    




def fast_are_synonyms(n, m, index_syn):

    try:
        index_syn[n]
    except:
        return False
    try:
        index_syn[m]
    except:
        return False
    return (index_syn[n] == index_syn[m])


def replace_in_sif(old_gene:str, new_gene:str, sif_content:str):
        sif_content = re.sub(r"^"+old_gene+"\t", r""+new_gene+"\t", sif_content, flags=re.MULTILINE)
        sif_content = re.sub(r"\t"+old_gene+"$", r"\t"+new_gene, sif_content, flags=re.MULTILINE)
        # sif_content = sif_content.replace('\t'+old_gene, '\t'+new_gene)
        return sif_content

def read_sif(sif: str) -> list:
    """
    Parses a SIF (Simple Interaction Format) string and returns a list of edges with signs.
    """
    edges = list()
    sif_content_list = sif.split("\n")
    # Remove last empty line
    sif_content_list = list(filter(None, sif_content_list))
    for edge in sif_content_list:
        source, sign, target = edge.split("\t")
        edges.append((source, target, {'sign': sign}))
    return edges

def identify_nodes(graph: nx.DiGraph) -> tuple:
    """
    Identifies input nodes, readout nodes, and complexes in the given graph.
    """
    inputs = set()
    readouts = set()

    for node in graph.nodes:
        if graph.in_degree(node) == 0:
            inputs.add(node)
        elif graph.out_degree(node) == 0:
            readouts.add(node)
    complexes = {out for _, out, sign in graph.edges.data('sign') if sign == 'PART_OF'}
    return inputs, readouts, complexes

def are_all_successors_in_list(node: str, successors_list: list, graph: nx.DiGraph) -> bool:
    """
    Checks if all successors of a node are present in the given list.
    """
    return all(successor in successors_list for successor in graph.successors(node))

def are_all_predecessors_in_list(node: str, predecessors_list: list, graph: nx.DiGraph) -> bool:
    """
    Checks if all predecessors of a node are present in the given list.
    """
    return all(predecessor in predecessors_list for predecessor in graph.predecessors(node))



def reduce_pkn(pkn_content: str, matrix_genes: list) -> tuple:
    """
    Reduces the given PKN content based on a list of matrix genes.
    """
    reduced_pkn_content = str()
    graph = nx.DiGraph(read_sif(pkn_content))
    
    to_keep = [n for n in graph.nodes if n in matrix_genes]

    complexes = list()
    # Remove complexes if their successors are not in the expression matrix
    for in_, out, sign in graph.edges.data('sign'):
        if sign == 'PART_OF':
            for succ in graph.successors(out):
                if succ in to_keep and out not in complexes:
                        complexes.append(out) 
 
    reduced_graph = graph.subgraph(to_keep+complexes)
    reduced_graph = reduced_graph.copy()
    inputs, readouts, _ = identify_nodes(reduced_graph)
    edges_to_remove = list()

    for in_, out, sign in reduced_graph.edges.data('sign'):
        if  (in_ in inputs and out in readouts) or \
            (in_ in inputs and out in complexes and are_all_successors_in_list(out, readouts, reduced_graph)) or \
            (in_ in complexes and out in readouts and are_all_predecessors_in_list(in_, inputs, reduced_graph)):
            edges_to_remove.append((in_, sign, out)) 
    for in_, _, out in edges_to_remove:
        reduced_graph.remove_edge(in_, out)

    for in_, out, sign in reduced_graph.edges.data('sign'):
        reduced_pkn_content += f'{in_}\t{sign}\t{out}\n'
    return reduced_pkn_content, edges_to_remove


def old_reduce_pkn(pkn_content:str, matrix_genes:list) -> str:
    reduced_pkn_content = str()
    graph = nx.DiGraph(read_sif(pkn_content))
    nodes = graph.nodes
    to_keep = [n for n in nodes if n in matrix_genes]
    complexes = list()
    for in_, out, sign in graph.edges.data('sign'):
        if sign == 'PART_OF':
            for succ in graph.successors(out):
                if succ in to_keep:
                    complexes.append(out) 
    reduced_graph = graph.subgraph(to_keep+complexes)
    reduced_graph_copy = reduced_graph.copy()
   
    inputs_with_one_succ = [node for node in reduced_graph if reduced_graph.in_degree(node) == 0 and reduced_graph.out_degree(node) == 1]
    readouts = [node for node in reduced_graph if reduced_graph.out_degree(node) == 0]
    deleted_inputs = list()
    for node in inputs_with_one_succ:
        succesor = next(reduced_graph.successors(node))
        if succesor in readouts: 
            deleted_inputs.append(node)
            reduced_graph_copy.remove_node(node)

    # reduced_graph = graph.subgraph(to_keep)
    for in_, out, sign in reduced_graph_copy.edges.data('sign'):
        # _,_,sign = reduced_graph_copy.edges[in_, out]['sign']
        reduced_pkn_content += f'{in_}\t{sign}\t{out}\n'
    return reduced_pkn_content, deleted_inputs




def reduce_expr_matrix(gene_expr_mtx_file: pd.DataFrame, genes_list:list, sif_content:str, annotation_len:int) -> set:
    """From the raw patrix, create a reduced matrix composed only of nodes 
        present also in the matrix and the PKN.

    Parameters
    ----------
    gene_expr_mtx_file : pd.DataFrame
        Raw expression matrix
    genes_list : list
        List of the PKN genes.
    sif_content : str
        SIF file content.

    Returns
    -------
    set
        Composed of the reduced matrix, modified SIF content file and some PKN genes informations compatible or not with genes in the matrix. 
    """
    modifications = {'pkn_name':[],
                    'matrix_name':[],
                    'new_pkn_name':[]}
    expr_mtx = load_data(gene_expr_mtx_file, index_name='Name')
    expr_mtx_genes = expr_mtx.columns
    to_keep = list()
    almost_in_expr_mtx = list()
    not_in_expr_mtx = list()
    perfect_matchs = list()
    modified_sif_content = sif_content

    for gene in genes_list:
        almost = False
        present = False
        for g in expr_mtx_genes:
            if gene == g: 
                to_keep.append(g)
                perfect_matchs.append(g)
                present=True
                break
            elif g.lower() == gene.lower(): 
                to_keep.append(g)
                modifications['pkn_name'].append(gene)
                modifications['matrix_name'].append(g)
                modifications['new_pkn_name'].append(gene+'__['+g+']')
                present = True
                modified_sif_content = replace_in_sif(gene, g, modified_sif_content)
                break
            else:
                if fast_are_synonyms(g, gene, index_syn):
                    to_keep.append(g)
                    modifications['pkn_name'].append(gene)
                    modifications['matrix_name'].append(g)
                    modifications['new_pkn_name'].append(gene+'__['+g+']')
                    present = True
                    modified_sif_content = replace_in_sif(gene, g, modified_sif_content)
                    break

                elif g.lower() in gene.lower():
                    almost = True
                

        if almost:
            almost_in_expr_mtx.append(gene)
        if  not almost and not present:
            not_in_expr_mtx.append(gene)


    annotation_columns = list(expr_mtx.columns[:annotation_len])
    to_keep = annotation_columns+to_keep
    reduced_expr_mtx = expr_mtx[to_keep]
    # to_keep = annotation_columns+perfect_matchs
    # perfect_matchs_reduced_mtx = expr_mtx[to_keep]
    # perfect_matchs_reduced_mtx = perfect_matchs_reduced_mtx.set_index('Name')
    return (reduced_expr_mtx, not_in_expr_mtx, almost_in_expr_mtx, modifications, modified_sif_content)


def plot(data, out_dir):
    

    Y_MIN_LIM = -0.5
    Y_MAX_LIM = 17

    X_MIN_LIM = 0
    X_MAX_LIM = 57

    COLOR_PLOT_STAGE = {
        'B1_B2':'#F5B3C4',
        'EPI':'#E9374E',
        'Morula':'#373C47',
        'PrE':'#85C258',
        'early_TE':'#65C3ED',
        'late_TE':'#B661C0',
        'medium_TE':'#5D6DAD',
        'ICM':'#F5EB6B'
    }



    pt_linregress = dict()

    empty_genes = list()

    stop = 0
    for gene in data.columns[5:]: # Skip the 5 firsts columns to have only genes
    # for gene in WGCNA_GS_LIST:
        # if stop==10:
        #     break
        # stop += 1




        # data = data.loc[(data[gene]> 0.5)]


        data_TE = data.loc[(data['clusterUmap']=='Morula') | (data['clusterUmap']=='B1_B2') | (data['clusterUmap']=='early_TE') | (data['clusterUmap']=='medium_TE') | (data['clusterUmap']=='late_TE')]
        data_PrE = data.loc[(data['clusterUmap']=='Morula') | (data['clusterUmap']=='B1_B2') | (data['clusterUmap']=='EPI') | (data['clusterUmap']=='PrE')]

        pt_linregress[gene] = dict()

        fig, (ax1, ax2) = plt.subplots(2, constrained_layout=True)

        TE_groups = data_TE.groupby('clusterUmap')
        PrE_groups = data_PrE.groupby('clusterUmap')
        
        if (len(data_TE[gene]) == 0) or (len(data_PrE[gene]) == 0): 
            empty_genes.append(gene)
            continue



        for name, group in TE_groups:
            ax2.plot(group.Pseudotime, group[gene], marker='o', linestyle='', color=COLOR_PLOT_STAGE[name], markersize=6, label=name)
        for name, group in PrE_groups:
            ax1.plot(group.Pseudotime, group[gene], marker='o', linestyle='', color=COLOR_PLOT_STAGE[name], markersize=6, label=name)
       
        linreg = stats.linregress(data_TE['Pseudotime'], data_TE[gene].astype(float))
        pt_linregress[gene]['TE'] = {
            'slope': linreg.slope,
            'intercept': linreg.intercept,
            'rvalue': linreg.rvalue,
            'pvalue': linreg.pvalue,
            'stderr': linreg.stderr,
            'intercept_stderr': linreg.intercept_stderr,
        }
        ax2.plot(data_TE['Pseudotime'], linreg.intercept + linreg.slope*data_TE['Pseudotime'], 'r', label='fitted line')


        linreg = stats.linregress(data_PrE['Pseudotime'], data_PrE[gene].astype(float))
        pt_linregress[gene]['PrE'] = {
            'slope': linreg.slope,
            'intercept': linreg.intercept,
            'rvalue': linreg.rvalue,
            'pvalue': linreg.pvalue,
            'stderr': linreg.stderr,
            'intercept_stderr': linreg.intercept_stderr,
        }
        ax1.plot(data_PrE['Pseudotime'], linreg.intercept + linreg.slope*data_PrE['Pseudotime'], 'r', label='fitted line')

       
        ax1.set_ylim(Y_MIN_LIM, Y_MAX_LIM)
        ax2.set_ylim(Y_MIN_LIM, Y_MAX_LIM)
        ax1.set_xlim(X_MIN_LIM, X_MAX_LIM)
        ax2.set_xlim(X_MIN_LIM, X_MAX_LIM)
        plt.suptitle(gene)
        ax1.legend( loc='center left', bbox_to_anchor=(1, 0.5), title='TE branch')
        ax2.legend( loc='center left', bbox_to_anchor=(1, 0.5), title='PrE branch')
        plt.savefig(out_dir+'plots/'+gene+'.pdf')
        plt.close()

    return (pt_linregress, empty_genes)



def analyze_linreg(data, output_name):

    

    output = {
        'reverse': list(),
        'diff_greater_than_same_sign' : {
            1.5: [],
            1: [],
            0.5: [],
            0.25: [],
            0.1: [],
            0.05: [],
            0.025: [],
            0.01: [],
        },
        'diff_greater_than_different_sign' : {
            1.5: [],
            1: [],
            0.5: [],
            0.25: [],
            0.1: [],
            0.05: [],
            0.025: [],
            0.01: [],
        }
    }

    values = [1.5, 1, 0.5, 0.25, 0.1, 0.05, 0.025, 0.01]

    for gene, g_data in data.items():
        pre = g_data['PrE']
        pre_slope = pre['slope']
        te = g_data['TE']
        te_slope = te['slope']

        pre_sign = '+' if pre_slope > 0 else '-' 
        te_sign = '+' if te_slope > 0 else '-' 

        if pre_sign != te_sign: 
            output['reverse'].append(gene)
            for val in values:
                if abs(pre_slope-te_slope) >= val :
                    output['diff_greater_than_different_sign'][val].append(gene)
        else:
            for val in values:
                if abs(pre_slope-te_slope) >= val :
                    output['diff_greater_than_same_sign'][val].append(gene)


    with open(output_name, 'w') as f:
        json.dump(output, f)


def get_input_genes_not_in_the_pkn(input, output):
    input = frozenset(input)
    output = frozenset(output)

    return list(input-output)


def remove_unknown_edges(sif_content:str) -> str:
    unknown_line_regex = '.*\sUNKNOWN\s.*\n'
    modified_sif_content = re.sub(unknown_line_regex, '', sif_content)
    return modified_sif_content



def run_pkn_analyze(sif_file:str, out_dir:str, gene_expr_mtx_file:str, input_genes_file:str, annotation_len:int):
    input_genes_list = read_file(input_genes_file)
    if out_dir[-1:] != '/': out_dir += '/'

    analyze_out_dir = f'{out_dir}analyze/'
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    if not os.path.exists(analyze_out_dir): os.makedirs(analyze_out_dir)
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    if not os.path.exists(analyze_out_dir+'plots/'): os.makedirs(analyze_out_dir+'plots/')

    out_linreg = 'analyze_lingreg.json'



    raw_data, genes_list = get_genes(sif_file)

    input_genes_not_in_the_pkn = get_input_genes_not_in_the_pkn(input_genes_list, genes_list)
    with open(sif_file, 'r') as f : 
        sif_content = "".join(f.readlines())
    reduced_expr_mtx, genes_not_in_expr_mtx, genes_almost_in_expr_mtx, modified_pkn_genes, modified_sif_content = reduce_expr_matrix(gene_expr_mtx_file, genes_list, sif_content, annotation_len)
    # Export reduced matrices
    reduced_expr_mtx.to_csv(out_dir+"pkn_gene_reduced_expr_mtx.csv")


    #remove Unknown edges
    modified_sif_content = remove_unknown_edges(modified_sif_content)

    # Reduce pkn
    reduced_pkn, removed_edges  = reduce_pkn(modified_sif_content, list(reduced_expr_mtx.columns))
    save_to_file(reduced_pkn, out_dir+'reduced_pkn.sif')
    # perfect_matchs_reduced_pkn, perfect_matchs_deleted_inputs = reduce_pkn(sif_content, list(perfect_matchs_reduced_mtx.columns))
    # save_to_file(perfect_matchs_reduced_pkn, out_dir+'perfect_matchs/reduced_pkn.sif')


    with open(analyze_out_dir+"input_genes_not_in_the_pkn.txt", "w") as outfile:
            for item in input_genes_not_in_the_pkn:
                outfile.write("{}\n".format(item))
    with open(analyze_out_dir+"genes_not_in_expr_mtx.txt", "w") as outfile:
            for item in genes_not_in_expr_mtx:
                outfile.write("{}\n".format(item))
    with open(analyze_out_dir+"genes_almost_in_expr_mtx.txt", "w") as outfile:
            for item in genes_almost_in_expr_mtx:
                outfile.write("{}\n".format(item))
    with open(analyze_out_dir+"removed_edges_from_the_pkn.txt", "w") as outfile:
                for item in removed_edges:
                    outfile.write("{}\n".format(item))
    # with open(analyze_out_dir+"perfect_matchs_deleted_inputs_from_the_pkn.txt", "w") as outfile:
    #             for item in perfect_matchs_deleted_inputs:
    #                 outfile.write("{}\n".format(item))
    with open(out_dir+'modified_pkn_genes.json', 'w') as outfile:
        json.dump(modified_pkn_genes, outfile, indent=4, sort_keys=True)

    

    # with open(out_dir+"modified_pkn.sif", 'w') as outfile:
    #     outfile.write(sif_content)

    # pt_linregress, empty_genes = plot(reduced_expr_mtx, analyze_out_dir)
    # with open(analyze_out_dir+"lin_regress.json", "w") as outfile:
    #         json.dump(pt_linregress, outfile, indent=4, sort_keys=True)
    # with open(analyze_out_dir+"empty_genes.txt", "w") as outfile:
    #         for item in empty_genes:
    #             outfile.write("{}\n".format(item))

    # analyze_linreg(pt_linregress, analyze_out_dir+out_linreg)



if __name__ == '__main__':


    
    sif_file = argv[1]
    out_dir = argv[2]
    gene_expr_mtx_file = argv[3]
    input_genes_file = argv[4]
    # no_predecessors_file = argv[5]
    # no_successors_file = argv[6]

    input_genes_list = read_file(input_genes_file)
    # no_predecessors = open(no_predecessors_file,'r').read().splitlines()
    # no_successors = open(no_successors_file,'r').read().splitlines()

    # EXP_NAME = 'Test_two_firsts_md3_fa'

    # sif_file = f'/home/e21g017n/Documents/work/Tools/pyBRAvo/MB_Test/{EXP_NAME}/-unified.sif'

    # gene_expr_mtx_file = '/home/e21g017n/Nextcloud/work/scRNASeq_data/data_analyze/data.zip'    

    if out_dir[-1:] != '/': out_dir += '/'

    analyze_out_dir = f'{out_dir}analyze/'
    if not os.path.exists(out_dir): os.makedirs(out_dir)
    if not os.path.exists(analyze_out_dir): os.makedirs(analyze_out_dir)
    if not os.path.exists(analyze_out_dir+'plots/'): os.makedirs(analyze_out_dir+'plots/')

    out_linreg = 'analyze_lingreg.json'


    index_std, index_syn = init_gene_synonyms_cache()

    raw_data, genes_list = get_genes(sif_file)

    input_genes_not_in_the_pkn = get_input_genes_not_in_the_pkn(input_genes_list, genes_list)
    with open(sif_file, 'r') as f : 
        sif_content = "".join(f.readlines())
    reduced_expr_mtx, genes_not_in_expr_mtx, genes_almost_in_expr_mtx, modified_pkn_genes, sif_content = reduce_expr_matrix(gene_expr_mtx_file, genes_list, sif_content)

    reduced_expr_mtx.to_csv(out_dir+"pkn_gene_reduced_expr_mtx.csv")

    with open(analyze_out_dir+"input_genes_not_in_the_pkn.txt", "w") as outfile:
            for item in input_genes_not_in_the_pkn:
                outfile.write("{}\n".format(item))
    with open(analyze_out_dir+"genes_not_in_expr_mtx.txt", "w") as outfile:
            for item in genes_not_in_expr_mtx:
                outfile.write("{}\n".format(item))
    with open(analyze_out_dir+"genes_almost_in_expr_mtx.txt", "w") as outfile:
            for item in genes_almost_in_expr_mtx:
                outfile.write("{}\n".format(item))

    with open(out_dir+'modified_pkn_genes.json', 'w') as outfile:
        json.dump(modified_pkn_genes, outfile, indent=4, sort_keys=True)

    with open(out_dir+"modified_pkn.sif", 'w') as outfile:
        outfile.write(sif_content)

    pt_linregress, empty_genes = plot(reduced_expr_mtx, analyze_out_dir)
    with open(analyze_out_dir+"lin_regress.json", "w") as outfile:
            json.dump(pt_linregress, outfile, indent=4, sort_keys=True)
    with open(analyze_out_dir+"empty_genes.txt", "w") as outfile:
            for item in empty_genes:
                outfile.write("{}\n".format(item))

    analyze_linreg(pt_linregress, analyze_out_dir+out_linreg)
