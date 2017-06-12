from ac_analysis.model.annotations import load_from_h5
from ac_analysis.lib.sequence_io.alignment_io import FastqExporter
import ac_analysis.model.ac_run_model as arm
import argparse

def filter_fastq_HQR(annotation_file):
    annotations = annotation_file
    output_name = (annotations + '_HQR.fastq')


    h5=load_from_h5(annotations)
    run_model=arm.ACRunModel()
    run_model.annotations=h5

    all_cells = h5.cells
    hqr,_ = h5.get_cell_annotation('high_quality_read')
    cells = all_cells[hqr]


    fq=FastqExporter()
    return fq.export_level_calls(run_model, output_name, cells=cells)


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, help='Input annotations file')
    args = parser.parse_args()

    return args



def main():
    args = get_args()


    filter_fastq_HQR(args.input)


if __name__ == '__main__':
    main()

