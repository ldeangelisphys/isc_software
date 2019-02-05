# ISC Analysis README
## Intro

This is a Python script to calculate Inter Subject Correlation starting from 4D nifti files from fMRI.

The correct running of the code requires python 3.X. To use on windows machines
we suggest installing anaconda. You also need to install the nibabel library, which is not 
by default in the anaconda distribution. To do so you can open the anaconda prompt and
type

> pip install nibabel


## Running the software
Once installed python, you can run the software by opening the anaconda prompt,
navigating to the location of the softare (cd /path/to/software) and by typing

> python calc_isc.py

In order to run correctly, the software needs a number of options to be
given by the user. To do so, we provide a settings file called settings.ini

## Setting File
In the file settings.ini the user can provide all the useful information about the input files that
need to be analyzed. These include the location of the input files.

> Folder with input data = C:/path/to/data/

This is the main folder containing all the data that will be needed for the analysis (all the Nifti files).

> Folder hierarchy = fld_S{XX}/preprocessed/ 

Inside the main folder, there could be a structure of subfolders, for example in order to determine the files
associated with each different subject to be analyzed. This structure can be specified as the Folder hierarchy,
where the symbol {XX} can be inserted whenever a number associated with the different subjects occurs.
The number of X here represent the number of digits used to represent the integer of each subject. For example,
for a structure like fld_S01/preprocessed/, fld_S02/preprocessed/, ..., we will use {XX}, whereas for
fld_S001/preprocessed/, fld_S002/preprocessed/, ..., we would use {XXX}. If all the files are in the folder with input
data specified before, this field can be left blank or set to '/'.

> N participants = 2

This is the number of subjects that we want to include in the analysis. Once set this number, the software
expects to find data starting from subject 1 to subject N participants.

> File name = 4D_niifile_S{XX}.nii

The name of the files with their extension should also be given to the software. Again, the symbol
{XX} can be used to replace the number associated to each subject, as in the folder hierarchy.

## Output

The software creates a folder where to store the generated output. This is created inside the folder with
input data given at the input and named ISC_XX-YY, where XX and YY are the first and last subject that
took part in the calculation. Inside this folder, the software saves one two .nii files for
each subject, one containing the leave-one-out ISC correlation and one its Z-Fisher transformation.
The software also stores the average signal of all the participants in this folder.
