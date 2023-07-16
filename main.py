from Bio import SeqIO
import pandas as pd
import numpy as np
import re
import os

MCSregionSequence = 'GAAACACCGACTTGCAGGTGCTTAAGGGATCCAGTATACTGGATCGATTGATCACACCTGCGGAT'
reference_file = pd.read_excel('Reference.xlsx')
matching_column = 'Sequence_Construct'
first_barcode_sequence = 'ACCG'
final_barcode_sequence = 'GTTT'
target_length = 20
barcode_score_length = 8
approximate_matching = False


def main():
    scans = pd.DataFrame()
    first_barcode_length = len(first_barcode_sequence)
    final_barcode_length = len(final_barcode_sequence)
    barcode_length_modifier_fr = barcode_score_length - first_barcode_length
    barcode_length_modifier_fi = barcode_score_length - final_barcode_length

    print(os.listdir('Sequencing data'))
    for files in os.listdir('Sequencing data'):
        if files.endswith('.ab1'):
            print(files)
            file_directory = os.path.join('Sequencing data', files)
            handle = open(file_directory, 'rb')
            record_iter = SeqIO.parse(handle, 'abi')
            for record in record_iter:
                well_location = record.annotations['sample_well'].decode()
                base_pair_list = list(record.seq)
                base_pair_list = pd.DataFrame(np.array(base_pair_list).reshape(len(base_pair_list), 1))
                sequence = record.seq
                phred_score = pd.DataFrame(record.letter_annotations)

                base_scores = pd.concat([base_pair_list, phred_score], axis=1)

                individual_bases = pd.DataFrame(base_pair_list, columns=['Individual Base'])
                phred_score['phred_percentage'] = phred_score['phred_quality'].apply(lambda x: 1 - (10 ** (-x / 10)))
                phred_score_mean = phred_score['phred_quality'].mean()
                final_dataframe = pd.concat([phred_score, individual_bases], axis=1)

                low_quality, med_quality, high_quality = final_dataframe[final_dataframe['phred_quality'] < 20]['phred_quality'], \
                                                         final_dataframe[(final_dataframe['phred_quality'] >= 20) & (
                                                                     final_dataframe['phred_quality'] < 40)]['phred_quality'], \
                                                         final_dataframe[final_dataframe['phred_quality'] >= 40]['phred_quality']

                low_quality_score, medium_quality_score, high_quality_score = (len(low_quality) / len(final_dataframe)) * 100, \
                                                                              (len(med_quality) / len(final_dataframe)) * 100, \
                                                                              (len(high_quality) / len(final_dataframe)) * 100

                low_quality_score, medium_quality_score, high_quality_score = format(low_quality_score, '.2f'), \
                                                                              format(medium_quality_score, '.2f'), \
                                                                              format(high_quality_score, '.2f')

                data = {'Sequence Barcode': [files],
                        'Well': [well_location],
                        'Sequence': [sequence],
                        'Sequence Quality Score': [phred_score_mean],
                        'HQuality %': [high_quality_score],
                        'MQuality %': [medium_quality_score],
                        'LQuality %': [low_quality_score],
                        'REGEX_Status': 'none',
                        'Beginning Quality Score': 0,
                        '20bp Sequence Quality Score': 0,
                        'Ending Quality Score': 0,
                        'Quality Score Ratio': 0
                        }

                sequence_regex = r'(%s(.{%s})%s)' % (str(first_barcode_sequence),
                                                     str(target_length),
                                                     str(final_barcode_sequence))

                mcs_regex = r'%s' % MCSregionSequence
                test = re.search(sequence_regex, str(sequence))
                mcs_test = re.search(mcs_regex, str(sequence))

                final_export = pd.DataFrame()
                if test is None:
                    if mcs_test is not None:
                        data['REGEX_Status'] = 'mcs'
                        data['Gene_Construct'] = 'MCS'
                        data['Sequence_Construct'] = 'MCS'
                        data['Matching status'] = 'MCS'
                        data = pd.DataFrame(data)
                        data.reset_index(inplace=True, drop=True)
                        final_export = data

                    else:
                        data['REGEX_Status'] = '-'
                        data = pd.DataFrame(data)
                        data.reset_index(inplace=True, drop=True)
                        final_export = data

                else:
                    matching_string = test.group(0)
                    matching_string = matching_string[first_barcode_length:-final_barcode_length]
                    data['REGEX_Status'] = matching_string

                    location_match = test.span()
                    transferred_scores = base_scores[location_match[0]
                                                     - barcode_length_modifier_fr:location_match[1]
                                                     + barcode_length_modifier_fi]
                    transferred_scores.reset_index(inplace=True, drop=True)

                    first_barcode = transferred_scores.loc[0:barcode_score_length - 1]

                    target_sequence = transferred_scores.loc[barcode_score_length:
                                                             (barcode_score_length + target_length) - 1]

                    second_barcode = transferred_scores.loc[barcode_score_length + target_length:
                                                            (barcode_score_length + target_length) +
                                                            (barcode_score_length-1)]

                    data['Beginning Quality Score'] = first_barcode['phred_quality'].mean()
                    data['20bp Sequence Quality Score'] = target_sequence['phred_quality'].mean()
                    data['Ending Quality Score'] = second_barcode['phred_quality'].mean()
                    barcode_quality_score = (data['Beginning Quality Score'] + data['Ending Quality Score']) / 2
                    data['Quality Score Ratio'] = data['20bp Sequence Quality Score'] / barcode_quality_score

                    data = pd.DataFrame(data).reset_index(drop=True)

                    if 0 < matching_string.count('N') <= 2:
                        if approximate_matching is True:
                            print(f'{matching_string} N detected in sequence')
                            query = matching_string.replace("N", ".")
                            print('finding closest match')

                            match_init = MatchClass(reference_file)
                            single_match_non_exact = match_init.reshape_dataframe(match_init.bool_match(query,
                                                                                                        matching_column))
                            if single_match_non_exact is not None:
                                final_export = pd.concat([data, single_match_non_exact], axis=1)
                                final_export['Matching status'] = 'match'

                    else:
                        match_init = MatchClass(reference_file)
                        single_match_exact = match_init.reshape_dataframe(match_init.bool_match(matching_string,
                                                                                                matching_column))
                        if single_match_exact is not None:
                            final_export = pd.concat([data, single_match_exact], axis=1)
                            final_export['Matching status'] = 'match'

                scans = pd.concat([scans, final_export])

    final_output = scans.replace(np.nan, 'no match', regex=True)
    text_file_output = final_output[['Sequence Barcode',
                                     'Well',
                                     'Matching status',
                                     'Gene_Construct',
                                     '20bp Sequence Quality Score']]
    text_file_output = text_file_output.rename(columns={'Sequence Barcode': 'Seq Barcode',
                                                        'Matching status': 'Status',
                                                        'Gene_Construct': 'Found',
                                                        '20bp Sequence Quality Score': 'Quality Score'})
    text_file_output['Manual Call'] = ''

    return text_file_output, final_output


class MatchClass:
    def __init__(self, matching_template):
        """
        initializes MatchClass
        :param matching_template: template to match query against
        :param matching_column: column within template to match query against
        """
        self.matching_template = matching_template

    def reshape_dataframe(self, *dataframe):
        """
        method that reshapes dataframe to cut out boolean score column and only return matching item row from template
        :param match_threshold: number that dictates scoring that qualifies as a match
        :param dataframe: variable that can hold multiple arguments that contain boolean dataframes
        :return: dataframe that has been reshape to only return matching row from template
        """

        def check_number_of_arguments(*args):
            """
            check number of arguments of *args variable in function instance
            :param args: Boolean dataframe
            :return: length of args list which shows how many arguments were given
            """
            return len(args)

        def convert_dataframe_from_bool(dataframe_input, true_as, false_as):
            """
            convert dataframe from bool into any user value; in this case, false is 0 and positive is 1
            :param dataframe_input: boolean dataframe to be inputted
            :param true_as: value to be used for true
            :param false_as: value to be used for false
            :return: converted dataframe
            """
            mask = dataframe_input.applymap(type) != bool
            d = {True: true_as, False: false_as}
            converted_dataframe = dataframe_input.where(mask, dataframe_input.replace(d)).T
            return converted_dataframe

        dataframe_check = check_number_of_arguments(*dataframe)
        if dataframe_check == 1:
            match_results = pd.DataFrame(dataframe).T
        else:
            match_results = pd.concat([*dataframe], axis=1)

        dataframe_slice = dataframe_check.__add__(1)

        match_results = convert_dataframe_from_bool(match_results, 1, 0)
        matching_sum = pd.DataFrame(np.sum(match_results, axis=0)).T
        match_results = pd.concat([match_results, matching_sum]).T

        final_output = pd.concat([match_results, self.matching_template], axis=1)

        database_export = final_output[(final_output[0] >= dataframe_check)]
        database_export = database_export.iloc[:, dataframe_slice:].reset_index(drop=True) # removes TRUE/FALSE dataframes and SUM column

        return database_export

    def bool_match(self, item_to_match, column_to_match):
        qr_string_match = self.matching_template.loc[:, column_to_match] == str(item_to_match)

        return qr_string_match


if __name__ == "__main__":
    cwd = os.getcwd()
    directory = cwd + r'\Sequencing data'

    text_file, final = main()

    text_file.to_csv("output.txt", sep="\t", index=None)
    final.to_excel("output.xlsx", index=None)
