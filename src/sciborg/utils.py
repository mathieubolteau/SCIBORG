#!/usr/bin/env python3
# -*- coding: utf-8 -*

try:
    import pandas as pd
    import csv
    import pkg_resources
    
except ImportError as E:
    print(E)
    print('Verify your dependancies.')
    exit()


def read_file(filename: str) -> list:
    """Read a file and return lines into a list.

    Parameters
    ----------
    filename : str
        Name of the file.

    Returns
    -------
    list
        List containing each line of the file.
    """
    with open(filename, 'r') as f:
        return f.read().splitlines()


def save_to_file(to_save: str, filename: str):
    """Save a string into file.

    Parameters
    ----------
    to_save : str
        String to save.
    filename : str
        Name of the file.
    """
    with open(filename, 'w') as f:
        f.write(to_save)


def list_to_file(list_to_save: list, filename: str):
    """Save a list into a file.

    Parameters
    ----------
    list_to_save : list
        List to save.
    filename : str
        Name of the file.
    """
    with open(filename, 'r') as f:
        for item in list_to_save:
            f.write(item)

def load_data(filename:str, index_name:str)-> pd.DataFrame:
    """Load a dataframe from a file.

    Parameters
    ----------
    filename : str
        Name of the file to import.

    Returns
    -------
    pd.DataFrame
        Loaded dataframe.
    """
    return pd.read_csv(filename,index_col=index_name)

def get_cell_class(cell:str, matrix:pd.DataFrame) -> str:
    return matrix.at[cell, 'clusterUmap']


def get_classes(matrix:pd.DataFrame) -> list:
    # TODO comments !!
    return matrix.clusterUmap.unique()

def init_gene_synonyms_cache():
    """

    :return:
    """
    index_syn = {}
    index_std = {}
    # with open(fullpath + '/Homo_sapiens.gene_info', newline='') as csvfile:
    synonyms_path = pkg_resources.resource_filename(__name__, 'pyBRAvo/src/bravo/Homo_sapiens.gene_info')
    
    # with open('/home/e21g017n/Nextcloud/work/gitlab_repos/pipeline/pipeline/pyBRAvo/src/bravo/Homo_sapiens.gene_info', newline='') as csvfile:
    with open(synonyms_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter='\t')
        next(reader)   # Skip first line
        for row in reader:
            index_std[row[2]] = row[4].split('|')

            index_syn[row[2]] = row[2]
            for syn in row[4].split('|'):
                index_syn[syn] = row[2]
    return index_std, index_syn


#print('--- Memory foot print cache ---')
index_std, index_syn = init_gene_synonyms_cache()