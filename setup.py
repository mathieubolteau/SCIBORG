from setuptools import setup, find_packages

# get version and other attributes
version_info = {}
with open("src/sciborg/version.py") as f:
    exec(f.read(), version_info)
with open("README.md", "r") as fh:
    long_description = fh.read()

# setup(
#     name = "SCIBORG"
#     description = "SCIBORG is an infering tool of Boolean models of stages involved in a cell differentiation system using singlel-cell transcriptomic data ."
#     readme = "README.md"
#     authors = [
#     { name = "Mathieu BOLTEAU", email = "mathieu.bolteau.pro@gmail.com" }
#     ]
#     entry_points={
#           'console_scripts': [
#               'SCIBORG=SCIBORG.SCIBORG:run',
#           ]
#       },
# )


# from setuptools import setup, find_packages

setup(
    name="sciborg",
    version=version_info['__version__'],
    author=version_info['__author__'],
    author_email=version_info['__author_email__'],
    description="SCIBORG is an infering tool of Boolean models of stages involved in a cell differentiation system using singlel-cell transcriptomic data ",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/mathieubolteau/SCIBORG",
    packages=find_packages(where='src',
                            include=['sciborg', 'sciborg.*', 'pyBRAvo', 'pyBRAvo.*']),
    package_dir={"":"src"},
    # packages=["sciborg"],
    # package_dir={"", "src"}
    # package_dir={"sciborg": "sciborg",
    #              "pyBRAvo":"sciborg/pyBRAvo",
    #              "tests":"sciborg/pyBRAvo/src/tests",
    #              "bravo":"sciborg/pyBRAvo/src/bravo",
    #              "pyBravo":"sciborg/pyBRAvo/src",},
    # packages=["sciborg", "pyBRAvo", "tests", "bravo", "pyBravo"],
    # install_requires=[
    #     "matplotlib>=3.7",
    #     "clingo>=5.5",
    #     "nxpd>=0.2",
    #     "networkx>=3.1",
    #     "SPARQLWrapper>=2.0",
    #     "requests>=2.31",
    #     "jupyter>=1.0",
    #     "rdflib>=6.3",
    #     "nose2>=0.14",
    #     "pandas>=2.0",
    #     "scipy>=1.7",
    #     "pydotplus>=2.0",
    #     "graphviz>=0.20",
    #     "clyngor>=0.4",
    # ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Operating System :: POSIX :: Linux",
        "Development Status :: 5 - Production/Stable",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    # packages = find_packages(),
    # package_data= {
    #     # all .dat files at any package depth
    #     '': ['**/*.gene_info'],
        
    # },
    entry_points={
          'console_scripts': [
              'sciborg=sciborg.sciborg:run',
          ]
      },
    include_package_data=True,
    #   package_data={
    #   'sciborg': ['**/*.gene_info'],
    #   'bravo': ['*.gene_info']
    #   }
      package_data={'sciborg':[
    "src/sciborg/**/*.lp",
    
    "src/sciborg/**/*.sh",
],
'pyBRAvo': ["src/bravo/Homo_sapiens.gene_info"]},
    data_files=[
        ('src/sciborg/pyBRAvo/src/bravo', ['src/sciborg/pyBRAvo/src/bravo/Homo_sapiens.gene_info']),
        ('src/sciborg/data', ['src/sciborg/data/pkn_construction/count_edges_type.lp', 
                              'src/sciborg/data/pkn_construction/get_complexes.lp',
                              'src/sciborg/data/pkn_construction/get_no_predecessors.lp',
                              'src/sciborg/data/pkn_construction/get_no_successors.lp',
                              'src/sciborg/data/pkn_construction/get_nodes.lp',
                              'src/sciborg/data/processing/get_ancestors.lp',
                              'src/sciborg/data/pseudo_perturbation_identification/problem.lp',
                              'src/sciborg/data/pseudo_perturbation_identification/run_pseudo_perturbation_identification.sh'
                              ]),
    ],
    #   package_data={'sciborg':}
    # data_files= [('', ["sciborg/pyBRAvo/src/bravo/Homo_sapiens.gene_info"])]
    
)