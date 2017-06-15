import os, sys
import pysam
from genomic_reports import summary
import matplotlib.pyplot as plt
import mpld3
from collections import defaultdict
import argparse


def read_fastq(fastq):
    seq_dict = {}
    f = open(fastq, 'r')
    for i, line in enumerate(f):
        if i%4 == 0:
            query_name = line.split('|')[1]
        if i%4 ==3:
            seq_length = len(line)
            seq_dict[query_name] = seq_length
    return seq_dict


def math_summary_individual_graph(path, amplicon_type, color_graph):
    amplicon_count_aggregate = {}
    readlenth_amplicons_aggregate = {}
    coverage_dict_aggregate = defaultdict(list)
    readlength_amplicons_fastq_aggregate = {}

    for filename in os.listdir(path):
        if filename.endswith('.fastq'):
            fastq = path + filename
            pokemon = filename.split('_')[4]
            bamfile = path + amplicon_type + '-' + pokemon + '_100_v5_MAPQ20.bam'
            coverage_txt = path + 'coverage_' + amplicon_type + '-' + pokemon + '_100_v5_MAPQ20.txt'
            bam_summary = summary.process_bam(bamfile)
            coverage_depth = open(coverage_txt, 'r')
            seq_lst = read_fastq(fastq)

            amplicon_count = {}
            readlenth_amplicons = {}
            coverage_dict = defaultdict(list)
            readlength_amplicons_fastq = {}

            for k,v in  bam_summary['reference_id'].iteritems():
                if v == '-1': continue
                if v not in amplicon_count_aggregate:
                    amplicon_count_aggregate[v] = 0
                amplicon_count_aggregate[v] +=1
                if v not in amplicon_count:
                    amplicon_count[v] = 0
                amplicon_count[v] += 1

            for line in coverage_depth:
                tab_split = line.split()
                amplicon = tab_split[0]
                base_num = tab_split[1]
                coverage = str(tab_split[2])
                coverage_dict_aggregate[amplicon].append(coverage)
                coverage_dict[amplicon].append(coverage)

            for k, v in coverage_dict.iteritems():
                total = 0
                for v1 in v:
                    total += int(v1)
                average = total / len(v)
                coverage_dict[k] = (average)


            for k, v in bam_summary['query_length'].iteritems():
                if bam_summary['reference_id'][k] != '-1':
                    if bam_summary['query_length'][k] not in readlenth_amplicons_aggregate:
                        readlenth_amplicons_aggregate[bam_summary['query_length'][k]] = 0
                    readlenth_amplicons_aggregate[bam_summary['query_length'][k]] += 1

                    if bam_summary['query_length'][k] not in readlenth_amplicons:
                        readlenth_amplicons[bam_summary['query_length'][k]] = 0
                    readlenth_amplicons[bam_summary['query_length'][k]] += 1

                    length_seq = seq_lst[bam_summary['query_name'][k].split('|')[1]]
                    if length_seq not in readlength_amplicons_fastq_aggregate:
                        readlength_amplicons_fastq_aggregate[length_seq] = 0
                    readlength_amplicons_fastq_aggregate[length_seq] += 1

                    if length_seq not in readlength_amplicons_fastq:
                        readlength_amplicons_fastq[length_seq] = 0
                    readlength_amplicons_fastq[length_seq] += 1



            fig = plt.figure()
            tPlot, axes = plt.subplots(
                nrows=4, ncols=1, sharex=False, sharey=False,
                gridspec_kw={'height_ratios': [10, 10, 10, 10]})

            ax1 = axes[0]
            ax2 = axes[1]
            ax3 = axes[2]
            ax4 = axes[3]

            ax1.bar(range(len(amplicon_count)), amplicon_count.values(), color=color_graph)
            ax1.set_xlabel('Amplicon Type')
            ax1.set_ylabel('Amplicon Count')
            ax1.set_title('Distribution of Amplicon')

            ax2.bar(range(len(coverage_dict)), coverage_dict.values(), color=color_graph)
            ax2.set_xlabel('Amplicon Type')
            ax2.set_ylabel('Depth of Coverage')
            ax2.set_title('Coverage Depth per Amplicon')

            lists2 = sorted(readlength_amplicons_fastq.items())  # sorted by key, return a list of tuples
            x2, y2 = zip(*lists2)  # unpack a list of pairs into two tuples
            ax3.plot(x2, y2, color=color_graph)
            ax3.set_xlabel('Read Length from Fastq file (bp)')
            ax3.set_ylabel('Amplicon Count')
            ax3.set_title('Distribution of Full Read Lengths')

            lists = sorted(readlenth_amplicons.items())  # sorted by key, return a list of tuples
            x, y = zip(*lists)  # unpack a list of pairs into two tuples
            ax4.plot(x, y, color=color_graph)
            ax4.set_xlabel('Read Length from BAM file(bp)')
            ax4.set_ylabel('Amplicon Count')
            ax4.set_title('Distribution of Sub-Read Lengths')

            tPlot.tight_layout()
            plots_html = mpld3.fig_to_html(tPlot)
            fig, ax = plt.subplots()

            output_name = bamfile[:-3] + '_analysis.html'
            with open((output_name), 'w') as output:
                 output.write(
                    '<b><font size ="5"> Amplicon Analysis Report: %s </b></font size><p>' % str(
                        'P-GA-' + pokemon + '_100_v5_MAPQ20'))
                 output.write(plots_html)

    for k, v in coverage_dict_aggregate.iteritems():
        total = 0
        for v1 in v:
            total += int(v1)
        average2 = total / len(v)
        coverage_dict_aggregate[k] = (average2)

    return amplicon_count_aggregate, readlenth_amplicons_aggregate, coverage_dict_aggregate, readlength_amplicons_fastq_aggregate


def plot_all_figs(path, amplicon_type, color_graph):

    amplicon_count, readlenth_amplicons, coverage_dict,readlength_amplicons_fastq = math_summary_individual_graph(path, amplicon_type, color_graph)
    fig = plt.figure()
    tPlot, axes = plt.subplots(
        nrows=4, ncols=1, sharex=False, sharey=False,
        gridspec_kw={'height_ratios': [7, 7, 7, 7] })

    ax1 = axes[0]
    ax2 = axes[1]
    ax3 = axes[2]
    ax4 = axes[3]


    ax1.bar(range(len(amplicon_count)), amplicon_count.values(), align="center", color = color_graph)
    ax1.set_xlabel('Amplicon Type')
    ax1.set_ylabel(' Amplicon Count')
    ax1.set_title('Distribution of Amplicon')


    ax2.bar(range(len(coverage_dict)), coverage_dict.values(), align="center", color = color_graph)
    ax2.set_xlabel('Amplicon Type')
    ax2.set_ylabel('Depth of Coverage')
    ax2.set_title('Average Coverage Depth per Amplicon')



    lists = sorted(readlength_amplicons_fastq.items())  # sorted by key, return a list of tuples
    x, y = zip(*lists)  # unpack a list of pairs into two tuples
    ax3.plot(x, y, color=color_graph)
    ax3.set_xlabel('Read Length from Fastq file (bp)')
    ax3.set_ylabel('Amplicon Count')
    ax3.set_title('Distribution of Full Read Lengths')

    lists2 = sorted(readlenth_amplicons.items())  # sorted by key, return a list of tuples
    x, y = zip(*lists2)  # unpack a list of pairs into two tuples
    ax4.plot(x, y, color=color_graph)
    ax4.set_xlabel('Read Length (bp)')
    ax4.set_ylabel('Amplicon Count')
    ax4.set_title('Distribution of Sub-Read Lengths')



    tPlot.tight_layout()
    plots_html = mpld3.fig_to_html(tPlot)
    fig, ax = plt.subplots()

    outputname = path + 'Analysis_Total_PETE_' + amplicon_type + '.html'
    with open((outputname), 'w') as output:
        output.write('<b><font size ="5"> Amplicon Analysis Report: Total PETE %s </b></font size><p>' % amplicon_type)
        output.write(plots_html)

    return


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='Location of input folder containing fastq, bam, txt.')
    parser.add_argument('--amp',type=str, help='type of amplicon that is in front of all files , e.g. P-GA' )
    parser.add_argument('--color', type=str, help='Color of graphs')

    args = parser.parse_args()

    return args


def main():
    args = get_args()

    plot_all_figs(args.input, args.amp, args.color)


if __name__ == '__main__':
    main()

#
# path = '/Users/burtone2/Documents/work/PETE_analysis/PETE_concat_data/PETE-SceI-GA/'
#
# print plot_all_figs(path,  'P-GA', 'red')