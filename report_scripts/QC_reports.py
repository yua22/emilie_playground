import os,sys
import matplotlib.pyplot as plt
import numpy as np
import mpld3
import itertools
from matplotlib import gridspec
import json
from pylab import *
from collections import defaultdict
from info_from_h5 import h5_parameter_table
import argparse


## getting expected reads from json file
def expected_read_seq(ExpState_json):
    json_file = open(ExpState_json, 'r')
    seq_attempt = json.load(json_file)[u'expData'][u'expectedReadSeq']
    return seq_attempt.strip('[').strip(']').strip('"').split(', "')

## values for key2
def report_calculations(path):

    ## groups is a #dictionary of dictionary
                #keys1: HQR, LQR
                #keys2: experiment number
                #values: apl_data(list), start_pos(list), num_reads, num_adapter_dimer, num_align_internal
    groups = defaultdict(dict)

    for filename in os.listdir(path):
        if filename.endswith('.csv'):
            cell_annotation_csv = filename
            experiment = filename[16:][:-3]
            ExpState_json = 'ExpState' + experiment + 'json'
            h5_file = path + 'annotations' + experiment + 'h5'

            expectedReadSeq = expected_read_seq(path + ExpState_json)
            cell_annotation_index = cell_annotation_csv.index('cell_annotations')
            experiment_num = cell_annotation_csv[cell_annotation_index + 17:][:-4]

            groups['HQR'][experiment_num] = [[], [], 0, 0, 0]
            groups['LQR'][experiment_num] = [[], [], 0, 0, 0]

            with open((path +cell_annotation_csv), 'r') as f:
                next(f)
                for line in f:
                    column = line.split(',')
                    align_procession_length = column[64]  # column BM:65, plus index = 64
                    high_quality_read = column[134]  # colulmn EE:135
                    align_template_index = column[69]  # column BR:70
                    align_edit_percent_identical = column[43]  # column AR:44
                    align_internal = column[57]  # BF: 58
                    align_aligned_template = column[29]  # AD: 30
                    fuctional_seq_pore = column[132]  # EC:133

                    if fuctional_seq_pore == 'True':
                        if high_quality_read == 'True':
                            groups['HQR'][experiment_num][2] +=1  #calculating number of HQR reads
                            groups['HQR'][experiment_num][0].append(int(align_procession_length)) #list of apl


                            try:
                                starting_positions = expectedReadSeq[int(align_template_index)].index(
                                    align_aligned_template.replace('-', '')[:len(expectedReadSeq[int(align_template_index)])])
                            except ValueError:
                                continue

                            groups['HQR'][experiment_num][1].append(starting_positions) #list starting lengths

                            if align_template_index == 0 and align_edit_percent_identical > 65:
                                groups['HQR'][experiment_num][3] += 1  #adapter_dimer
                            if align_internal == 'True':
                                groups['HQR'][experiment_num][4] += 1 # internal_alignment




                        else:
                            groups['LQR'][experiment_num][2] += 1  # calculating number of HQR reads
                            groups['LQR'][experiment_num][0].append(int(align_procession_length))  # list of apl
                            try:
                                starting_positions = expectedReadSeq[int(align_template_index)].index(
                                    align_aligned_template.replace('-', '')[
                                    :len(expectedReadSeq[int(align_template_index)])])
                            except ValueError:
                                continue

                            groups['LQR'][experiment_num][1].append(starting_positions)  # list starting lengths

                            if align_template_index == 0 and align_edit_percent_identical > 65:
                                groups['LQR'][experiment_num][3] += 1  # adapter_dimer
                            if align_internal == 'True':
                                groups['LQR'][experiment_num][4] += 1  # internal_alignment

    return groups, h5_file

def single_report(path, output_path):

    lst_single_report = []
    groups, h5_file = report_calculations(path)
    for k,v in groups.iteritems():
        for k1,v1 in v.iteritems():
            if k1 not in lst_single_report:
                lst_single_report.append(k1)

    for i in lst_single_report:

        ## writing plots to html format
        fig = plt.figure()

        gs1 = gridspec.GridSpec(2, 1)
        ax1 = fig.add_subplot(gs1[0])
        ax2 = fig.add_subplot(gs1[1])

        ## plots

        if groups['HQR'][i][2] != 0:
            binwidth = 50  # chosen by eye
            ax1.hist([groups['HQR'][i][0], groups['LQR'][i][0]], bins=binwidth, stacked=True, color=['grey', 'blue'],
                     label=['HQR', 'LQR'])
            ax1.legend()
            ax1.set_xlabel('Align Procession Length')
            ax1.set_ylabel('Frequency')
            ax1.set_title('Histogram of Align Procession Length for ' + str(i))
        else:
            binwidth = 50  # chosen by eye
            ax1.hist([groups['LQR'][i][0]], bins=binwidth, stacked=True, label=['LQR'])
            ax1.legend()
            ax1.set_xlabel('Align Procession Length')
            ax1.set_ylabel('Frequency')
            ax1.set_title('Histogram of Align Procession Length for ' + str(i))

        sorted_start_pos_data_HQR = np.sort(groups['HQR'][i][1])
        sorted_start_pos_data_LQR = np.sort(groups['LQR'][i][1])

        yvals = np.arange(len(sorted_start_pos_data_HQR)) / float(len(sorted_start_pos_data_HQR)-1)
        yvals2 = np.arange(len(sorted_start_pos_data_LQR)) / float(len(sorted_start_pos_data_LQR)-1)
        ax2.plot(sorted_start_pos_data_HQR, yvals, color='grey', label='HQR')
        ax2.plot(sorted_start_pos_data_LQR, yvals2, color='blue', label='LQR')
        ax2.set_ylim([0.95, 1])
        ax2.set_xlabel('Start Positions')
        ax2.set_ylabel('Cumulative Frequency')
        ax2.legend(loc='lower right')
        ax2.set_title('Cumulative Frequency of Start Positions Summary')

        gs1.tight_layout(fig)

        plots_html = mpld3.fig_to_html(fig)

        fig, ax = plt.subplots()

        if groups['HQR'][i][2] != 0:
            percent_total_adapter_dimer_HQR = float(
                "{0:.4f}".format(groups['HQR'][i][3] / float(groups['HQR'][i][2]))) * 100
            percent_total_align_internal_HQR = float(
                "{0:.4f}".format(groups['HQR'][i][4] / float(groups['HQR'][i][2]))) * 100

        else:
            percent_total_adapter_dimer_HQR = 0.0
            percent_total_align_internal_HQR = 0.0


        if groups['LQR'][i][2] == 0 :
            percent_total_adapter_dimer_LQR = 0.0
            percent_total_align_internal_LQR = 0.0

        else:
            percent_total_adapter_dimer_LQR = float(
                "{0:.4f}".format(groups['LQR'][i][3] / float(groups['LQR'][i][2]))) * 100
            percent_total_align_internal_LQR = float(
                "{0:.4f}".format(groups['LQR'][i][4] / float(groups['LQR'][i][2]))) * 100


        total_adapter_dimer_HQR = str(str(groups['HQR'][i][3]) + ' (' + str(percent_total_adapter_dimer_HQR) + '%)')
        total_align_internal_HQR = str(str(groups['HQR'][i][4]) + ' (' + str(percent_total_align_internal_HQR) + '%)')
        total_adapter_dimer_LQR = str(str(groups['LQR'][i][3]) + ' (' + str(percent_total_adapter_dimer_LQR) + '%)')
        total_adapter_internal_LQR = str(str(groups['LQR'][i][4]) + ' (' + str(percent_total_align_internal_LQR) + '%)')



        tbl = [[' ', 'HQR', 'LQR'],
               ['num of adapter dimers:', total_adapter_dimer_HQR, total_adapter_dimer_LQR],
               ['num of internal alignments:', total_align_internal_HQR, total_adapter_internal_LQR]]

        ## code for table to html
        cols = ["<td>{0}</td>".format("</td><td>".join(t)) for t in tbl]

        # then use it to join the rows (tr)
        rows = "<tr>{0}</tr>".format("</tr>\n<tr>".join(cols))

        with open((output_path + str(i) +'.html'), 'w') as output:
            output.write(('<b><font size ="5"> Analysis Report: %s  </b></font size><p>' % str(i)))
            output.write(h5_parameter_table(h5_file, 'Experiment Parameters'))
            output.write("""<HTML>
                                         <h1><font size ="3">HQR vs LQR list </font size></h1>
                                          <table>
                                            {0}
                                          </table><br>

                                </HTML>""".format(rows)
                         )

            output.write(plots_html)

def summary_HQR_data(path, output_path):


    fig = plt.figure()
    gs1 = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs1[0])
    ax2 = fig.add_subplot(gs1[1])

    legend_label = []
    groups, h5_file = report_calculations(path)


    total_HQR_align_procession = []
    label_HQR_align_procession = []
    store_tables_parameters = []
    store_table_adapters= []


    for k,v in groups['HQR'].iteritems():
        total_HQR_align_procession.append(v[0])
        label_HQR_align_procession.append(k)

        legend_label.append(k)

        sorted_start_pos_data_HQR = np.sort(v[1])
        yvals = np.arange(len(sorted_start_pos_data_HQR)) / float(len(sorted_start_pos_data_HQR)-1)
        ax2.plot(sorted_start_pos_data_HQR, yvals)
        ax2.legend(legend_label)
        ax2.set_ylim([0.97, 1])
        ax2.set_xlabel('Start Positions')
        ax2.set_ylabel('Cumulative Frequency')
        ax2.legend(loc='lower right')
        ax2.set_title('Cumulative Frequency of Start Positions Summary HQR')

        ## parameters
        h5_table = h5_parameter_table(h5_file, ('Experiment Parameters: ' + k))
        store_tables_parameters.append(h5_table)



        if v[2] == 0: continue
        percent_total_adapter_dimer_HQR = float(
            "{0:.4f}".format(v[3] / float(v[2]))) * 100
        percent_total_align_internal_HQR = float(
            "{0:.4f}".format(v[4] / float(v[2]))) * 100

        total_adapter_dimer_HQR = str(str(v[3]) + ' (' + str(percent_total_adapter_dimer_HQR) + '%)')
        total_align_internal_HQR = str(str(v[4]) + ' (' + str(percent_total_align_internal_HQR) + '%)')


        tbl = [[' ', 'HQR'],
               ['num of adapter dimers:', total_adapter_dimer_HQR],
               ['num of internal alignments:', total_align_internal_HQR]]

        ## code for table to html
        cols = ["<td>{0}</td>".format("</td><td>".join(t)) for t in tbl]

        # then use it to join the rows (tr)
        rows = "<tr>{0}</tr>".format("</tr>\n<tr>".join(cols))

        table_HQR = ("""<HTML>
                            <h1><font size ="3">HQR values  </font size></h1>
                                   <table>
                                        {0}
                                   </table><br>

                    </HTML>""".format(rows)
                )

        store_table_adapters.append(table_HQR)



    ax1.hist(total_HQR_align_procession, bins=50, alpha=0.7, stacked = True, label = label_HQR_align_procession)
    ax1.legend()
    ax1.set_xlabel('Align Procession Length')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Histogram of HQR Align Procession Length')



    gs1.tight_layout(fig)
    plots_html = mpld3.fig_to_html(fig)
    fig, ax = plt.subplots()



    with open((output_path + 'Summary_HQR_per_experiment.html'), 'w') as output:
        output.write('<b><font size ="5"> Analysis Report: Summary of HQR data </b></font size><p>')
        for p,a in zip(store_tables_parameters, store_table_adapters):
            output.write(p)
            output.write(a)
        output.write(plots_html)


    return

def summary_all_data(path, output_path):
    groups, h5_file = report_calculations(path)

    total_apl_HQR = []
    total_start_pos_HQR = []
    total_num_read_HQR = 0
    total_num_adapter_dimer_HQR = 0
    total_num_align_internal_HQR = 0

    total_apl_LQR = []
    total_start_pos_LQR = []
    total_num_read_LQR = 0
    total_num_adapter_dimer_LQR = 0
    total_num_align_internal_LQR = 0




    for k,v in groups['HQR'].iteritems():
            total_apl_HQR.extend(v[0])
            total_start_pos_HQR.extend(v[1])
            total_num_read_HQR += v[2]
            total_num_adapter_dimer_HQR += v[3]
            total_num_align_internal_HQR += v[4]

    for k, v in groups['LQR'].iteritems():
        total_apl_LQR.extend(v[0])
        total_start_pos_LQR.extend(v[1])
        total_num_read_LQR += v[2]
        total_num_adapter_dimer_LQR += v[3]
        total_num_align_internal_LQR += v[4]



    ## writing plots to html format
    fig = plt.figure()

    gs1 = gridspec.GridSpec(2, 1)
    ax1 = fig.add_subplot(gs1[0])
    ax2 = fig.add_subplot(gs1[1])


    binwidth = 50  # chosen by eye
    ax1.hist([total_apl_HQR, total_apl_LQR], bins=binwidth, stacked=True, color = ['grey', 'blue'], label=['HQR', 'LQR'])
    ax1.legend()
    ax1.set_xlabel('Align Procession Length')
    ax1.set_ylabel('Frequency')
    ax1.set_title('Histogram of Align Procession Length Summary' )

    sorted_start_pos_data_HQR = np.sort(total_start_pos_HQR)
    sorted_start_pos_data_LQR = np.sort(total_start_pos_LQR)

    yvals = np.arange(len(sorted_start_pos_data_HQR)) / float(len(sorted_start_pos_data_HQR))
    yvals2 = np.arange(len(sorted_start_pos_data_LQR)) / float(len(sorted_start_pos_data_LQR))
    ax2.plot(sorted_start_pos_data_HQR, yvals, color = 'grey', label = 'HQR')
    ax2.plot(sorted_start_pos_data_LQR, yvals2, color = 'blue', label = 'LQR')
    ax2.set_ylim([0.95, 1])
    ax2.set_xlabel('Start Positions')
    ax2.set_ylabel('Cumulative Frequency')
    ax2.legend(loc ='lower right')
    ax2.set_title('Cumulative Frequency of Start Positions Summary')

    gs1.tight_layout(fig)

    plots_html = mpld3.fig_to_html(fig)

    fig, ax = plt.subplots()

    ## tables

    if total_num_read_HQR == 0:
        percent_total_adapter_dimer_HQR = 0.0
        percent_total_align_internal_HQR = 0.0

    else:
        percent_total_adapter_dimer_HQR = float("{0:.4f}".format(total_num_adapter_dimer_HQR / float(total_num_read_HQR))) * 100
        percent_total_align_internal_HQR = float("{0:.4f}".format(total_num_align_internal_HQR / float(total_num_read_HQR))) * 100

    if total_num_read_LQR == 0:
        percent_total_adapter_dimer_LQR = 0.0
        percent_total_align_internal_LQR = 0.0
    else:
        percent_total_adapter_dimer_LQR = float("{0:.4f}".format(total_num_adapter_dimer_LQR / float(total_num_read_LQR))) * 100
        percent_total_align_internal_LQR = float("{0:.4f}".format(total_num_align_internal_LQR / float(total_num_read_LQR))) * 100


    HQR_adapter_dimers = str((str(total_num_adapter_dimer_HQR)+ ' ('+ str(percent_total_adapter_dimer_HQR) + '%)'))
    LQR_adapter_dimers = str((str(total_num_adapter_dimer_LQR) + ' (' + str(percent_total_adapter_dimer_LQR) + '%'))
    HQR_internal_alignments = str((str(total_num_align_internal_HQR) + ' (' + str(percent_total_align_internal_HQR) + '%)'))
    LQR_internal_alignments = str((str(total_num_align_internal_LQR)+ ' (' + str(percent_total_align_internal_LQR )+ '%)'))

    tbl = [[' ', 'HQR' , 'LQR',],
           ['num of adapter dimers: ', HQR_adapter_dimers, LQR_adapter_dimers],
           ['num of internal alignments:', HQR_internal_alignments, LQR_internal_alignments ]]

    ## code for table to html
    cols = ["<td>{0}</td>".format("</td><td>".join(t)) for t in tbl]

    # then use it to join the rows (tr)
    rows = "<tr>{0}</tr>".format("</tr>\n<tr>".join(cols))

    ## code for table to html
    cols = ["<td>{0}</td>".format("</td><td>".join(t)) for t in tbl]

    # then use it to join the rows (tr)
    rows = "<tr>{0}</tr>".format("</tr>\n<tr>".join(cols))

    with open((output_path + 'Summary.html'), 'w') as output:
        output.write('<b><font size ="5"> Analysis Report: Aggregate Summary </b></font size><p>')

        output.write("""<HTML>
                                      <h1><font size ="3">HQR vs LQR list </font size></h1>
                                       <table>
                                         {0}
                                       </table><br>

                             </HTML>""".format(rows)
                     )

        output.write(plots_html)

    return output

def get_all(path, output_path):
    print single_report(path, output_path)
    print summary_HQR_data(path,output_path)
    print summary_all_data(path, output_path)

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='Location of input folder containing csv, json, h5 files.')
    parser.add_argument('--output', type=str, help='Location of output folder ')
    args = parser.parse_args()

    return args


def main():
    args = get_args()

    get_all(args.input, args.output)


if __name__ == '__main__':
    main()


