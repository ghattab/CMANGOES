import os
import contextlib
import pandas as pd
import multiprocessing as mp
import cmangoes
import tqdm
from timeit import default_timer as timer


path_datasets = os.path.join('..', 'Data', 'Original_datasets')
list_of_datasets = [
    'ace_vaxinpad',
    'acp_anticp',
    'acp_iacp',
    'acp_mlacp',
    'afp_amppred',
    'afp_antifp',
    'aip_aippred',
    'aip_antiinflam',
    'amp_antibp',
    'amp_antibp2',
    'amp_csamp',
    'amp_fernandes',
    'amp_gonzales',
    'amp_iamp2l',
    'amp_modlamp',
    'atb_antitbp',
    'atb_iantitb',
    'avp_amppred',
    'avp_avppred',
    'bce_ibce',
    'cpp_cellppd',
    'cpp_cellppdmod',
    'cpp_cppredfl',
    'cpp_kelmcpp',
    'cpp_mixed',
    'cpp_mlcpp',
    'cpp_mlcppue',
    'cpp_sanders',
    'hem_hemopi'
]


def encode_one_dataset_parallel(path_one_dataset):
    path_one_dataset_fasta = os.path.join(path_one_dataset, 'seqs.fasta')
    path_one_dataset_classes = os.path.join(path_one_dataset, 'classes.txt')
    path_one_dataset_output_dir = os.path.join(path_one_dataset, 'csv', 'all')
    os.makedirs(path_one_dataset_output_dir, exist_ok=True)
    path_one_dataset_smiles = os.path.join(
                path_one_dataset_output_dir, 'resulting_smiles.smi')

    binary_encoding_values = [True, False]
    center_encoding_values = [True, False]
    level = 12
    smiles_file_exists_flag = False
    smiles_list = []

    for binary_encoding in binary_encoding_values:
        for center_encoding in center_encoding_values:
            output_distinct_name = 'binary_' if binary_encoding\
                else 'discretized_'
            output_distinct_name += 'centered_' if center_encoding\
                else 'shifted_'
            output_distinct_name += 'levels_1_and_2_'

            # If we don't have SMILES file, we create it using FASTA
            if not smiles_file_exists_flag:
                with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                    smiles_list = cmangoes.convert_fasta_to_smiles(
                        path_one_dataset_fasta, path_one_dataset_smiles)
                num_of_lines = len(smiles_list)
                names = range(1, num_of_lines + 1)
                smiles_file_exists_flag = True

            path_one_data_set_output_encoding = os.path.join(
                path_one_dataset_output_dir,
                output_distinct_name + 'encoding.csv')

            with open(os.devnull, "w") as f, contextlib.redirect_stdout(f):
                finalized_encoding = cmangoes.encode_molecules(
                    smiles_list, names, binary_encoding=binary_encoding,
                    center_encoding=center_encoding, level=level,
                    output_path=path_one_dataset_output_dir)

                cmangoes.csv_export(
                    finalized_encoding,
                    output_path=path_one_data_set_output_encoding)

            # We open the encodings file and append one column at the end that
            # represents classes, labeled 'y'. We use the classes.txt file
            # given at the data dir. Additionaly, we add one column at the
            # beginning with sequence identifiers 'Seq_#', without any label of
            # the column.

            # Creating the last column with classes
            one_dataset_classes = []
            with open(path_one_dataset_classes, 'r') as classes_file:
                one_dataset_classes = [
                    int(one_line.rstrip()) for one_line in classes_file]

            one_dataset_num_of_encodings = len(one_dataset_classes)

            # Creating the first column with sequence labels
            one_dataset_first_column_labels = [
                'Seq_' + str(i) for i in range(
                    1, one_dataset_num_of_encodings+1)]

            # Appending both columns to the exported CMANGOES CSV file with
            # encodings
            one_dataset_encodings_df = pd.read_csv(
                path_one_data_set_output_encoding)

            # Appending the first column
            one_dataset_encodings_df.insert(
                0, '', one_dataset_first_column_labels)

            # Appending the last column
            one_dataset_encodings_df['y'] = one_dataset_classes

            # Saving the dataframe on the same location, with the same name
            one_dataset_encodings_df.to_csv(
                path_one_data_set_output_encoding, index=False)

            # We release the memory the dataframe holds
            tmp_lst = [one_dataset_encodings_df]
            del one_dataset_encodings_df
            del tmp_lst

    # We remove the SMILES file that was created for the encoding of one
    # dataset by cmangoes
    if os.path.exists(path_one_dataset_smiles):
        os.unlink(path_one_dataset_smiles)

    return None


def run_parallel():
    # TODO: Generalize run_sequential as well. Use tqdm for progress
    # instead of the print function
    num_of_datasets = len(list_of_datasets)
    print('==================================================================')
    print('\t\t\tPARALLEL')
    print('Encoding ' + str(num_of_datasets) + ' datasets using CMANGOES')
    print('==================================================================')

    NUMBER_OF_PROCESSES = mp.cpu_count() - 2
    path_all_datasets = [os.path.join(path_datasets, one_dataset)
                         for one_dataset in list_of_datasets]

    with mp.Pool(processes=NUMBER_OF_PROCESSES) as p:
        for _ in tqdm.tqdm(
                p.imap_unordered(encode_one_dataset_parallel,
                                 path_all_datasets),
                total=num_of_datasets):
            pass

    print('==================================================================')
    print('\t\t\tPARALLEL ENCODING IS DONE')
    print('==================================================================')

    return None


def run_sequential(bool_flag_time):
    # Traversing all datasets and encoding them using CMANGOES
    num_of_datasets = len(list_of_datasets)
    dataset_counter = 1
    print('==================================================================')
    print('Encoding ' + str(num_of_datasets) + ' datasets using CMANGOES')
    print('==================================================================')

    if bool_flag_time:
        list_time_results = []

    for one_dataset in list_of_datasets:
        print('Progress ' + str(dataset_counter) + '/' + str(num_of_datasets),
              '|', one_dataset)
        path_one_dataset = os.path.join(path_datasets, one_dataset)
        path_one_dataset_fasta = os.path.join(path_one_dataset, 'seqs.fasta')
        path_one_dataset_classes = os.path.join(path_one_dataset,
                                                'classes.txt')
        path_one_dataset_output_dir = os.path.join(path_one_dataset, 'csv',
                                                   'all')
        os.makedirs(path_one_dataset_output_dir, exist_ok=True)
        path_one_dataset_smiles = os.path.join(
                    path_one_dataset_output_dir, 'resulting_smiles.smi')

        binary_encoding_values = [True, False]
        center_encoding_values = [True, False]
        level = 12
        smiles_file_exists_flag = False
        smiles_list = []

        for binary_encoding in binary_encoding_values:
            for center_encoding in center_encoding_values:
                output_distinct_name = 'binary_' if binary_encoding\
                    else 'discretized_'
                output_distinct_name += 'centered_' if center_encoding\
                    else 'shifted_'
                output_distinct_name += 'levels_1_and_2_'

                print('\nCurrent encoding:', one_dataset, output_distinct_name)
                print('======================================================')

                # If we don't have SMILES file, we create it using FASTA
                if not smiles_file_exists_flag:
                    smiles_list = cmangoes.convert_fasta_to_smiles(
                        path_one_dataset_fasta, path_one_dataset_smiles)
                    num_of_lines = len(smiles_list)
                    names = range(1, num_of_lines + 1)
                    smiles_file_exists_flag = True

                path_one_data_set_output_encoding = os.path.join(
                    path_one_dataset_output_dir,
                    output_distinct_name + 'encoding.csv')

                if bool_flag_time:
                    timer_start = timer()

                finalized_encoding = cmangoes.encode_molecules(
                    smiles_list, names, binary_encoding=binary_encoding,
                    center_encoding=center_encoding, level=level,
                    output_path=path_one_dataset_output_dir)

                if bool_flag_time:
                    timer_end = timer()
                    list_time_results.append(timer_end - timer_start)

                cmangoes.csv_export(
                    finalized_encoding,
                    output_path=path_one_data_set_output_encoding)

                # We open the encodings file and append one column at the end
                # that represents classes, labeled 'y'. We use the classes.txt
                # file given at the data dir. Additionaly, we add one column at
                # the beginning with sequence identifiers 'Seq_#', without any
                # label of the column.

                # Creating the last column with classes
                one_dataset_classes = []
                with open(path_one_dataset_classes, 'r') as classes_file:
                    one_dataset_classes = [
                        int(one_line.rstrip()) for one_line in classes_file]

                one_dataset_num_of_encodings = len(one_dataset_classes)

                # Creating the first column with sequence labels
                one_dataset_first_column_labels = [
                    'Seq_' + str(i) for i in range(
                        1, one_dataset_num_of_encodings+1)]

                # Appending both columns to the exported CMANGOES CSV file with
                # encodings
                one_dataset_encodings_df = pd.read_csv(
                    path_one_data_set_output_encoding)

                # Appending the first column
                one_dataset_encodings_df.insert(
                    0, '', one_dataset_first_column_labels)

                # Appending the last column
                one_dataset_encodings_df['y'] = one_dataset_classes

                # Saving the dataframe on the same location, with the same name
                one_dataset_encodings_df.to_csv(
                    path_one_data_set_output_encoding, index=False)

                # We release the memory the dataframe holds
                tmp_lst = [one_dataset_encodings_df]
                del one_dataset_encodings_df
                del tmp_lst

        # We remove the SMILES file that was created for the encoding of one
        # dataset by cmangoes
        if os.path.exists(path_one_dataset_smiles):
            os.unlink(path_one_dataset_smiles)

        dataset_counter += 1
        print('==============================================================')

    return list_time_results


def test_performance(int_number_of_runs):
    bool_flag_time = True
    path_experiment_data = os.path.join('..', 'Data',
                                        'Performance_experiments')

    list_encoding_names = []
    binary_encoding_values = [True, False]
    center_encoding_values = [True, False]

    for string_one_dataset in list_of_datasets:
        for binary_encoding in binary_encoding_values:
            for center_encoding in center_encoding_values:
                output_distinct_name = '_binary_' if binary_encoding\
                    else '_discretized_'
                output_distinct_name += 'centered_' if center_encoding\
                    else 'shifted_'
                output_distinct_name += 'levels_1_and_2_'

                list_encoding_names.append(
                    string_one_dataset + output_distinct_name)

    df_results = pd.DataFrame()
    df_results['Encodings'] = list_encoding_names

    for i in int_number_of_runs:
        list_results_of_run = run_sequential(bool_flag_time)
        df_results['Run_' + str(i)] = list_results_of_run
        df_results.to_csv(os.path.join(path_experiment_data, 'results.csv'),
                          index=False)

    return None


def main():
    # run_parallel()

    int_number_of_runs = 3
    test_performance(int_number_of_runs)

    return None


if __name__ == '__main__':
    main()
