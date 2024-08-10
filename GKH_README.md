# Demultiplexing
This repo takes a series of Illumina reads barcoded using dual-matched barcodes, and filters them based on the quality score (>20), barcode identity (present in the file of barcodes used), and whether they are matched. 

Results for the data in `/projects/bgmp/shared/2017_sequencing/` in tsv format is in this repo, under [this](results_summarized/) folder. Using this quality score cutoff, the data was split as follows:

| Category       |  Frequency  |
| -------------- | :---------: |
| Matched count: | 330,738,415 |
|  Hopped count: |   679,459   |
| Unknown count: | 31,828,861  |
These results are also summarized in this figure:
[fig1](results_summarized/fig1.png)

[Assignment the first](Assignment-the-first/) comprises the outline of the algorithm and the main functions, as well as exploratory data analysis. Commands and notes can be found [here](Assignment-the-first/Assignment_The_First_Working_Notes.pdf). Quality score cutoff was determined by plotting the data using a script that can be run easily with [run.sh](run.sh). This was used to generate the following figure:
[fig0](Assignment-the-first/quality_hist.py)

[Assignment the third](Assignment-the-third/) contains the actual code, and appropriate tests.
Unit testing of individual functions is done by [this](Assignment-the-third/test_functions.sh) script, and functional testing is done by [this](Assignment-the-third/test_whole_script.sh). Test files are included in the repo. 

Versions can be found [here](Assignment-the-third/Versions.md).

