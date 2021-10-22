import os
import pathlib
import argparse
import math
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib import colors
from pysmiles import read_smiles
from IPython.display import clear_output
from openbabel import openbabel
from Bio import SeqIO


def convert_fasta_to_smiles(input_fasta_file_path, output_smiles_file_path):
    obConversion = openbabel.OBConversion()
    obConversion.SetInAndOutFormats('fasta', 'smi')

    mol = openbabel.OBMol()

    with open(input_fasta_file_path, 'r') as input_file,\
         open(output_smiles_file_path, 'a') as output_file:
        for i, fasta_string in enumerate(SeqIO.parse(input_file, 'fasta')):
            obConversion.ReadString(mol, str(fasta_string.seq))
            output_smiles_string = obConversion.WriteString(mol)
            # print(i+1, '-', output_smiles_string)
            output_file.write(output_smiles_string)

    print('Successfully converted FASTA into SMILES.')

    return None


def create_graph_for_molecule(mol):
    adj = nx.adjacency_matrix(mol, weight='order').todense()

    return nx.from_numpy_matrix(adj)


def get_labels_from_elements(elements):
    labels = {}
    for idx, el in elements:
        labels[idx] = "{}: {}".format(idx, el)

    return labels


# This function was used for testing purposes
def plot_molecule_graph(G, labels):
    pos = nx.spring_layout(G)
    nx.draw(G, pos=pos, node_size=400)
    nx.draw_networkx_labels(G, pos, labels, font_size=10)
    plt.show()

    return None


def encode_molecule(mol, plot_molecule=False):
    elements = mol.nodes(data="element")
    G = create_graph_for_molecule(mol)

    if (plot_molecule):
        labels = get_labels_from_elements(elements)
        plot_molecule_graph(G, labels)

    carbons = [c for c in elements if c[1].lower() == "c"]
    neighborhoods = {}

    for idx, c in carbons:
        neighbors_idx = list(G[idx].keys())
        neighbors_elements = [elements[key] for key in neighbors_idx]
        neighborhoods["level{}".format(idx)] = pd.Series(neighbors_elements)

    return pd.DataFrame.from_dict(neighborhoods)


def next_perfect_square(N):
    nextN = math.floor(math.sqrt(N)) + 1

    return nextN * nextN


def center_matrix(m, target_dim):
    cur_dim = list(m.shape)
    steps = target_dim - cur_dim[0]
    for i in range(steps):
        if i % 2 == 0:
            to_append = np.zeros((cur_dim[1], 1))
            m = np.append(to_append, m, axis=1)
            cur_dim[0] += 1
            to_append = np.zeros((1, cur_dim[0]))
            m = np.append(to_append, m, axis=0)
            cur_dim[1] += 1
        else:
            to_append = np.zeros((1, cur_dim[0]))
            m = np.append(m, to_append, axis=0)
            cur_dim[1] += 1
            to_append = np.zeros((cur_dim[1], 1))
            m = np.append(m, to_append, axis=1)
            cur_dim[0] += 1

    return(m)


def shift_matrix(m, target_dim):
    cur_dim = list(m.shape)
    steps = target_dim - cur_dim[0]
    for i in range(steps):
        to_append = np.zeros((1, cur_dim[0]))
        m = np.append(m, to_append, axis=0)
        cur_dim[1] += 1
        to_append = np.zeros((cur_dim[1], 1))
        m = np.append(m, to_append, axis=1)
        cur_dim[0] += 1

    return(m)


def get_unique_atoms(mol):
    atoms = mol.nodes(data="element")
    unique_atoms = set()
    for atom_tuple in atoms:
        unique_atoms.add(atom_tuple[1])

    return unique_atoms


# Function to generate dummy encoding of smiles strings
def dummy_encode_molecules(smiles, binary_encoding=True, print_progress=False,
                           plot_molecule=False):
    res = []
    number_of_elements = len(smiles)

    if not binary_encoding:
        unique_atoms = set()

    if print_progress:
        progress = 0

    for molecule in smiles:
        if print_progress:
            clear_output(wait=True)
            progress += 1
            print("encoding molecule {} of {}".format(
                progress, number_of_elements))

        mol = read_smiles(molecule, explicit_hydrogen=True)

        if not binary_encoding:
            unique_atoms.update(get_unique_atoms(mol))

        encoding = encode_molecule(mol, plot_molecule=plot_molecule)
        dummy_encoding = pd.get_dummies(encoding)

        if not binary_encoding:
            c_columns = [col for col in dummy_encoding if col.endswith("C")]
            n_columns = [col for col in dummy_encoding if col.endswith("N")]
            o_columns = [col for col in dummy_encoding if col.endswith("O")]
            p_columns = [col for col in dummy_encoding if col.endswith("P")]
            s_columns = [col for col in dummy_encoding if col.endswith("S")]

            dummy_encoding[c_columns] = dummy_encoding[c_columns] * 2
            dummy_encoding[n_columns] = dummy_encoding[n_columns] * 3
            dummy_encoding[o_columns] = dummy_encoding[o_columns] * 4
            dummy_encoding[p_columns] = dummy_encoding[p_columns] * 5
            dummy_encoding[s_columns] = dummy_encoding[s_columns] * 6

        res.append(dummy_encoding)

    return res


# Function to normalize dummy encoding
def normalize_encodings(dummy_encodings, names, center_encoding=True):

    max_dim = 0
    squared_matrices = []
    output_dict = {}

    for dummy in dummy_encodings:
        dummies_as_list = dummy.transpose().values.tolist()
        dummies_flat = [item for sublist in dummies_as_list
                        for item in sublist]
        filler_list = [0] * (
            next_perfect_square(len(dummies_flat)) - len(dummies_flat))
        encoding_squared = dummies_flat + filler_list
        dimension = int(math.sqrt(len(encoding_squared)))
        max_dim = max(max_dim, dimension)
        squared_matrices.append(
            np.array(encoding_squared).reshape(dimension, dimension))

    if center_encoding:
        print('''Maximum dimension is {}x{}, centering smaller matrices in
                 {}x{} grid...'''.format(max_dim, max_dim, max_dim, max_dim))
    else:
        print('''Maximum dimension is {}x{}, shifting smaller encoding to match
                 maximum dimension...'''.format(max_dim, max_dim))

    for i in range(len(dummy_encodings)):
        if center_encoding:
            centered = center_matrix(squared_matrices[i], max_dim)
            output_dict[names[i]] = list(np.ravel(centered.astype(int)))

        else:
            shifted = shift_matrix(squared_matrices[i], max_dim)
            output_dict[names[i]] = list(np.ravel(shifted.astype(int)))

    return output_dict


# Function to generate images from normalized encoding
def generate_imgs_from_encoding(normalized_encoding, binary_encoding=True,
                                folder_name="encoding_images",
                                print_progress=False):

    if print_progress:
        clear_output(wait=True)
        progress = 0
        number_of_items = len(normalized_encoding)

    for name, encoding in normalized_encoding.items():

        if print_progress:
            clear_output(wait=True)
            progress += 1
            print("generating image {} of {}".format(
                progress, number_of_items))

        plt.figure(figsize=(10, 10))
        # plt.title(name, fontsize=26)
        ax = plt.gca()
        ax.axes.xaxis.set_ticks([])
        ax.axes.yaxis.set_ticks([])

        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_color("grey")
            ax.spines[axis].set_linewidth(3)

        dim = int(math.sqrt(len(encoding)))

        # Path seperator was originally presented as '\'
        # FIX: os.path.join is used instead, for uniformity
        # dirname = os.path.realpath(".") + ("\{}".format(folder_name))

        dirname = os.path.join(os.path.realpath("."), folder_name)
        if not os.path.exists(dirname):
            os.makedirs(dirname)
        # filename = dirname + ("\{}.png".format(name))
        filename = os.path.join(dirname, name)

        if binary_encoding:
            cmap = colors.ListedColormap(["lightgrey", "black"])
            cmap_bounds = [0, 0.1, 1]
            norm = colors.BoundaryNorm(cmap_bounds, cmap.N)

        else:
            cmap = colors.ListedColormap(["lightgrey", "white", "black",
                                          "blue", "red", "orange", "yellow"])
            cmap_bounds = [0, 1, 2, 3, 4, 5, 6, 7]
            norm = colors.BoundaryNorm(cmap_bounds, cmap.N)

        plt.imshow(np.array(encoding).reshape(dim, dim), cmap=cmap, norm=norm)
        plt.savefig(filename)
        plt.close()

    # print("Saved images to folder ./{}".format(folder_name))
    print("Saved images to folder ." + os.path.sep() + folder_name)

    return None


# Wrapper function for dummy, normalize, image generation (optional)
def encode_molecules(
        smiles, names, binary_encoding=True, center_encoding=True,
        plot_molecule=False, print_progress=False, generate_images=False,
        output_path="encoding_images"):

    dummies = dummy_encode_molecules(
        smiles=smiles, binary_encoding=binary_encoding,
        print_progress=print_progress, plot_molecule=plot_molecule)
    normalized_encoding = normalize_encodings(
        dummy_encodings=dummies, names=names, center_encoding=center_encoding)
    if generate_images:
        generate_imgs_from_encoding(
            normalized_encoding=normalized_encoding,
            binary_encoding=binary_encoding, folder_name=output_path,
            print_progress=print_progress)

    return normalized_encoding


# CSV export of normalized encoding
def csv_export(normalized_encoding, classes=pd.DataFrame(),
               output_path="encoding.csv"):
    encoding_as_df = pd.DataFrame.from_dict(
        normalized_encoding, orient="index")
    encoding_as_df = encoding_as_df.reset_index(drop=True)

    if not len(classes.index) == 0:
        classes = classes.reset_index(drop=True)
        encoding_as_df = encoding_as_df.join(classes)

    encoding_as_df.to_csv(output_path, index=False)
    print("Exported to ", output_path)

    return None


# Generate encodings and export CSVs
# Helper function to generate all permutatations of encodings
def generate_all_encodings(smiles, names, data_set_identifier,
                           classes=pd.DataFrame()):

    # Paths were hard-coded before. Below is the proper definition
    binary_centered_out_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_binary_centered_imgs")
    binary_centered_csv_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_binary_centered.csv")
    binary_shifted_out_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_binary_shifted_imgs")
    binary_shifted_csv_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_binary_shifted.csv")
    discretized_centered_out_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_discretized_centered_imgs")
    discretized_centered_csv_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_discretized_centered.csv")
    discretized_shifted_out_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_discretized_shifted_imgs")
    discretized_shifted_csv_path = os.path.join(
        "..", "Images", data_set_identifier, data_set_identifier,
        "_discretized_shifted.csv")

    print("Generating binary centered encoding...")
    binary_centered = encode_molecules(
        smiles, names, print_progress=True, generate_images=True,
        output_path=binary_centered_out_path)
    csv_export(
        binary_centered, classes=classes, output_path=binary_centered_csv_path)

    print("Generating binary shifted encoding...")
    binary_shifted = encode_molecules(
        smiles, names, print_progress=True, center_encoding=False,
        generate_images=True, output_path=binary_shifted_out_path)
    csv_export(
        binary_shifted, classes=classes, output_path=binary_shifted_csv_path)

    print("Generating discretized centered encoding...")
    discretized_centered = encode_molecules(
        smiles, names, binary_encoding=False, print_progress=True,
        generate_images=True,
        output_path=discretized_centered_out_path)
    csv_export(
        discretized_centered, classes=classes,
        output_path=discretized_centered_csv_path)

    print("Generating discretized shifted encoding...")
    discretized_shifted = encode_molecules(
        smiles, names, binary_encoding=False, center_encoding=False,
        print_progress=True, generate_images=True,
        output_path=discretized_shifted_out_path)
    csv_export(
        discretized_shifted, classes=classes,
        output_path=discretized_shifted_csv_path)

    return None


# This function was once a part of the main function
# It is used to create sample data sets that are mentioned in the paper
def create_datasets():
    # Paths were hard-coded before. Below is the proper definition
    amino_acid_path = os.path.join("..", "Data", "amino_acids",
                                   "amino_acids.csv")
    ace_vaxinpad_path = os.path.join("..", "Data", "ace_vaxinpad",
                                     "ace_vaxinpad.smiles")
    ace_vaxinpad_classes_path = os.path.join("..", "Data", "ace_vaxinpad",
                                             "ace_vaxinpad_classes.txt")
    hiv_protease_path = os.path.join("..", "Data", "hiv_protease",
                                     "hiv_protease.smiles")
    hiv_protease_classes_path = os.path.join("..", "Data", "hiv_protease",
                                             "hiv_protease_classes.txt")

    # Read and prepare data
    amino_acids = pd.read_csv(amino_acid_path, delimiter=",", header=None,
                              names=["names", "smiles"])
    amino_acids_smiles = list(amino_acids["smiles"])
    amino_acids_names = list(amino_acids["names"])

    ace_vaxinpad = pd.read_csv(ace_vaxinpad_path, delimiter=",", header=None,
                               names=["smiles"])
    ace_vaxinpad_smiles = list(ace_vaxinpad["smiles"])
    ace_vaxinpad_names = list(range(1, len(ace_vaxinpad_smiles) + 1))
    ace_vaxinpad_classes = pd.read_fwf(ace_vaxinpad_classes_path, header=None,
                                       names=["y"])

    hiv_protease = pd.read_csv(hiv_protease_path, delimiter=",", header=None,
                               names=["smiles"])
    hiv_protease_smiles = list(hiv_protease["smiles"])
    hiv_protease_names = list(range(1, len(hiv_protease_smiles) + 1))
    hiv_protease_classes = pd.read_fwf(hiv_protease_classes_path, header=None,
                                       names=["y"])

    # Amino Acids
    generate_all_encodings(
        smiles=amino_acids_smiles, names=amino_acids_names,
        data_set_identifier="amino_acids")

    # Ace vaxinpad
    generate_all_encodings(
        smiles=ace_vaxinpad_smiles, names=ace_vaxinpad_names,
        data_set_identifier="ace_vaxinpad", classes=ace_vaxinpad_classes)

    # HIV Protease
    generate_all_encodings(
        smiles=hiv_protease_smiles, names=hiv_protease_names,
        data_set_identifier="hiv_protease", classes=hiv_protease_classes)

    return None


def main():

    program_name = 'cmangoes'
    program_description = '''cmangoes: Carbon-based Multi-level Atomic
                             Neighborhood Encodings'''

    input_help = 'A required path-like argument.'
    encoding_help = 'A required character argument.'
    level_help = 'An optional integer argument.'
    image_help = 'An optional boolean argument. Default: False'
    padding_help = 'An optional character argument. Default: c'
    graph_help = 'An optional integer argument.'

    output_default = os.path.join('.', 'CMANGOES_Results')
    output_help = '''An optional path-like argument. For parsed paths, the
                     directory must exist beforehand.
                     Default: ''' + output_default

    input_error = 'Input file path is bad or the file does not exist'
    graph_error = '''Graph should be an integer >=1 and <=number of sequences
                     in the input file'''
    output_error = '''Output directory path is bad or the directory does not
                      exist'''

    argument_parser = argparse.ArgumentParser(
        prog=program_name, description=program_description)

    # Adding arguments
    allowed_encodings = ['b', 'd']
    allowed_paddings = ['c', 's']
    allowed_levels = [1, 2, 3]

    argument_parser.add_argument('input_file', type=pathlib.Path,
                                 help=input_help)
    argument_parser.add_argument('encoding', type=str, help=encoding_help,
                                 choices=allowed_encodings)
    argument_parser.add_argument('--level', type=int, help=level_help,
                                 choices=allowed_levels)
    argument_parser.add_argument('--image', default=False, help=image_help)
    argument_parser.add_argument('--padding', type=str, default='c',
                                 help=padding_help, choices=allowed_paddings)
    argument_parser.add_argument('--show_graph', type=int, help=graph_help)
    argument_parser.add_argument('--output_dir', type=pathlib.Path,
                                 help=output_help)

    # Parsing arguments
    arguments = argument_parser.parse_args()

    # Additional argument inspection
    if not os.path.exists(arguments.input_file):
        argument_parser.error(input_error)

    # TODO: Implement the check for <= num_of_sequences in the input file
    if arguments.show_graph is not None:
        if arguments.show_graph <= 0:
            argument_parser.error(graph_error)

    if arguments.output_dir is not None:
        if not os.path.exists(arguments.output_dir):
            argument_parser.error(output_error)

    return None


# Input
# - fasta
# - smiles

# Argument
# - level 1 :3
# - encoding b: binary and d: discretized
# - image 1,0
# - padding c: centered and s: shifted
# - user path for output write
# - show-graph id : show intermediary graph for molecule index

# Output

# - smiles as output written to file (iff fasta used as input) SMI or smiles
# - binary or discretized files CSV
# - images filename as index from encoding CSV
# - log 


if __name__ == '__main__':
    main()
