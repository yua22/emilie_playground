README REPORTS 

Reports:
QC_reports.py
Amplicon_reports.py
________________________________
QC_reports.py  

Goal:
To generate intuitive and clean reports of QC metrics to determine which templates and parameters works well for Sample Prep team.  

Graphs and Tables:
Tables: - experiment parameters (token, mode, notes, questions)
	  - HQR vs LQR 
Graphs: - histogram of align procession length
	   - cumulative frequency of start positions 


Conclusion:
HTML reports are created with interactive graphs (from mpld3) and links to cumulus to view parameters. They can be generated using a terminal command: 
python QC_reports.py --input /path_to_where_csv_json_h5_files_are/ --output /path_to_desired_location_of_output_files/

Methods: 
Inputs necessary: 
•	‘ExpState_run_name.json’
•	‘cell_annotations_run_name.csv’
•	‘annotations_run_name.h5’

Note, script will not run if run_name is not identical. Ie:
ExpState_170413_SAM_01_zapdos_WAU17R11C15.json, 
cell_annotations_170413_SAM_01_zapdos_WAU17R11C15.csv, 
annotations_170413_SAM_01_zapdos_WAU17R11C15.h5 


Requirements: info_from_h5.py 


 Output/Results : 
	•Html report for each individual run   —  170413_SAM_01_zapdos_WAU17R11C15.html
	•Html report summarizing all HQRs      — Summary_HQR_per_experiment.html 
	•Html report aggregating all runs in folder  — Summary.html 


_________________________________ 

Amplicon_reports.py


Goal:
To generate intuitive and clean reports for Amplicon analysis for Sample Prep team.  

Graphs: - Distribution of Amplicons 
	- Coverage Depth per Amplicon 
	- Distribution of Full Read Lengths
	- Distribution of Sub-Read Lengths 

Conclusion:
HTML reports are created with interactive graphs(mpld3). They can be generated using a terminal command:  
python Amplicon_reports.py --input /path_to_where_csv_json_h5_files_are/ --output /path_to_desired_location_of_output_files/


Methods: 
Inputs necessary: 
•	‘homo_date_SAM_run_name.annotations.h5_HQR.fastq’
•	‘amplicon_type_run_name_v5_MAPQ20.bam’
•	‘coverage_amplicon_tupe_run_name_v5_MAPQ20.txt’

Note, script will not run if run_name is not identical. Ie. 
homo_170525_SAM_01_goldeen_WAL12R11C07.annotations.h5_HQR.fastq, 
P-SceI-goldeen_100_v5_MAPQ20.bam, 
coverage_P-SceI-goldeen_100_v5_MAPQ20.txt

Txt file is generated using samtools ‘samtools depth’. Parameters used: samtools depth -a -d 500000 


 Output/Results : 
	•Html report for each individual run - Analysis_P-SceI-goldeen_100_v5_MAPQ20.html 
	•Html report aggregating all runs in folder  — Analysis_Total_PETE-SceI-Y.html








