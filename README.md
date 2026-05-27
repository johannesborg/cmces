Compile cmces.cpp with gcc or g++. It is necessary to link the boost graph library (https://www.boost.org/doc/libs/1_75_0/libs/graph/doc/index.html) and the openbabel (library https://open-babel.readthedocs.io/en/latest/index.html).

Input for the program is a csv file where each line contains comma separated SMILES strings for which you wish to find the maximum common subgraphs.

Output from the program is a csv file, where each line is a comma separated list of SMILES strings, which are the maximum common subgraphs for each corresponding line in the input file.

The program is run as follows:

./cmces inputfile.csv outputfile.csv option1 option2

Here option1 one must be either "prune" or "no-prune", where "prune" activates the extra pruning of type-0 edges in the modular product, and "no-prune" uses the classical product.

option2 must be either "order" or "no-order", where order activates the ordering heuristic with the minmax similarity and "no-order" uses the ordering provided in the input file.

The script order_graphs.py is used for ordering input graphs with graph kernel methods.

The input for the script is a csv file where each line contains comma separated SMILES strings

The output is a csv file where each line contains the same SMILES strings as in the input, just ordered differently.

The script is run as follows:

python order_graphs.py inputfile outputfile kernel

Where kernel is an integer in [1,2,3], with 1: vertexhistogram kernel, 2: NSPD kernel, 3: WL-OA kernel.
