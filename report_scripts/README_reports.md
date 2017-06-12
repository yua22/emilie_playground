README REPORTS 

QC_reports.py  

Goal:
To generate intuitive and clean reports of QC metrics to determine which templates and parameters works well for Sample Prep team.  

Conclusion:
An HTML report is created with interactive graphs and links to cumulus to view parameters. They can be generated using a terminal command:  ‘
python reports_v3.py --input /path_to_where_csvjson_h5_files_are/ --output /path_to_where_file_for_output_files_should_be_  /

Methods: 
Inputs necessary: 
•	‘ExpState_run_name.json’
•	‘cell_annotations_run_name.csv’
•	‘annotations_run_name.h5’

Note, script will not run if run_name is not identical. Ie. ExpState_170413_SAM_01_zapdos_WAU17R11C15.json, 
cell_annotations_170413_SAM_01_zapdos_WAU17R11C15.csv 
annotations_170413_SAM_01_zapdos_WAU17R11C15.h5 


Requirements: info_from_h5.py 


 Output/Results : 
	•	Html report for each individual run.  —  170413_SAM_01_zapdos_WAU17R11C15.html
	•	Html report summarizing all HQRs   — Summary_HQR_per_experiment.html 
	•	Html report aggregating all runs in folder  — Summary.html 








