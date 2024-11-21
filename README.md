# SCIBORG - using Single-Cell data to Infer BOolean networks modeling Regulation of Genes

SCIBORG is an inferring tool of Boolean models of stages involved in a cell differentiation system using single-cell transcriptomic data 


The framework searches to compute families of Boolean networks that are both compatible with scRNA-seq data and prior regulatory knowledge. It is composed of three steps:
1. The PKN reconstruction.
2. The experimental design construction.
3. The BN inference.

### Prerequisites

SCIBORG only works under Python 3.8. Two submodules, included in the SCIBORG package, are required:
* [pyBRAvo](https://github.com/mathieubolteau/pyBRAvo) 
* [caspo](https://github.com/mathieubolteau/caspo)

### Installation

1. Clone the repository including submodules.
    ```sh 
    git clone --recursive git@github.com:mathieubolteau/SCIBORG.git
    ```
2. Install anaconda/miniconda if not already done to create conda env.
   Please follow instructions here: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
3. Create a virtual env, **under Python 3.8**.
    ```sh 
    cd SCIBORG
    conda env create --file environment.yml
    ```
4. Activate your virtual env.
   ```sh 
   conda activate sciborg_env
   ```
5. Install clyngor package (an automatic solution is under development) 
   ```sh 
   pip install clyngor
   ```
6. Install SCIBORG package
    ```sh 
    pip install .
    ```
7. Verify that SCIBORG is installed
   ```sh 
   sciborg --help
   ```

## Usage

```sh
usage: sciborg [-h] [-pkn] [-pp] [-po] [-bn] configuration_file

positional arguments:
  configuration_file    configuration file path

optional arguments:
  -h, --help            show this help message and exit
  -pkn, --pkn-construction
                        run the PKN construction
  -pp, --pseudo-perturbation-identification
                        run the pseudo-perturbation identification
  -po, --pseudo-observation-diff-maxi
                        run the pseudo-observation difference maximization
  -bn, --bn-inference   run the BN inference
```

SCIBORG works with a configuration file to fix the parameters of each step of the framework. 

To make the *PKN reconstruction step* run the following command: 
```
sciborg --pkn-construction <CONFIG FILE>
```

To make the *Experimental design construction* run the following command. Notice that step required that the previous step has already been completed.
```
sciborg --pseudo-perturbation-identification --pseudo-observation-diff-maxi <CONFIG_FILE>
```
Note that the sub-steps can be run independently.

To make the *BN inference* step, run the following command:
```
sciborg --bn-inference <CONFIG_FILE>
```

## Demo

In the `demo` folder, a toy dataset to test the installation of SCIBORG and see the expected results for this toy dataset.

To reproduce the expected results, using the configuration file present in the `demo/data` folder, run the following steps.

1. Modify *endpoint* link in the `demo_config.ini` file. See *PKN CONSTRUCTION* section.

2. Run the following command:

```
sciborg --pkn-construction --pseudo-perturbation-identification --pseudo-observation-diff-maxi --bn-inference ./demo/data/demo_config.ini
```

The execution should take 3/4 minutes and produce results similar to those in `demo/expected_results` folder.

## Author
* Mathieu Bolteau 
    * [@mathieubolteau](https://github.com/mathieubolteau)
    * [mathieu.bolteau.pro@gmail.com](mailto:mathieu.bolteau.pro@gmail.com)
    * Nantes Université, École Centrale Nantes, CNRS, LS2N, UMR 6004, F-44000, Nantes