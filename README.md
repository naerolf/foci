# Foci clustering tools
Summary: clustering for foci dStorm

## To check the interactive notebook

Just open the IPython Notebook _foci_data_ :wink:

## For image processing

Call the python script _foci_process.py_ with the following syntaxis:
```
python foci_process.py <input_file> <output_folder> <minimum_cluster_size> <cluster_probability_threshold> 
```
Example:
```
python .\foci_process.py 'D:\RSI\icastro\foci\data\28.1_BLEO_1000_1.txt.csv' 'D:\RSI\icastro\foci\output' 300 0.5
```