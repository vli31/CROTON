<a id='sec0'></a>
# CROTON
Please follow the following instructions to host CROTON on your local computer. CROTON is an automated and variant-aware CRISPR/Cas9 editing outcome predictor built on deep multi-task convolutional neural networks and neural architecture search. 

CROTON predicts the following CRISPR/Cas9 editing outcome statistics: (1) 1 bp insertion frequency, (2) 1 bp deletion frequency, (3) deletion frequency, (4) 1 bp frameshift frequency, (5) 2 bp frameshift frequency, and (6) total frameshift frequency. 

See https://doi.org/10.1093/bioinformatics/btab268 for more information.

### Get the latest code
Clone the Github Repository. If you have a previous version, make sure to pull the latest commits:
```
git clone https://github.com/vli31/CROTON
cd CROTON
git pull
```
If you see `Already up to date` in the terminal, it means the code you have is the latest version.

## Installation

CROTON is available as a [Django](https://www.djangoproject.com/) webserver, in which users input a 60 bp CRISPR/Cas9 target sequence into a query box. Note that the cut site should be in the middle of the input sequence, with 30 nucleotides to the left and right of the double stranded break.

CROTON can be installed with the following steps:

### 1. Creating an Anaconda microenvironment
To run CROTON, it is easiest to use a conda microenvironment. If you do not have Anaconda, you can download it at this link: https://www.anaconda.com/products/individual. 

The environment used for CROTON is available in this repository in the file 'croton.yml'. It is recommended that you use the same versions of Django (3.0.3), Keras (2.2.4) and TensorFlow (1.14.0).
```
conda env create -f croton.yml
conda activate croton
```
Now, you should be in the 'croton' conda microenvironment. Make sure you see `(croton)` instead of `(base)` to the left of the '$' in your terminal window.

### 2. Setting up Django 
To set up Django make sure you are in the directory 'CROTON' and have activated the microenvironment 'croton,' then run:
```
python manage.py runserver
```
Now, the development server should be available at http://127.0.0.1:8000/. Click on the link to access CROTON.

You can quit the server with CONTROL-C.

## Testing

The CROTON web interface should look like this:

<img pull-left src="croton_screenshot.png">

If you input the example nucleotide sequence "TCCAGGGCCTAATCTGACCGTCCTAGATACCTCAGGGTGGGCAATACGAGGTAATGGCAG," the webserver should output:

- 1 bp Insertion Probability: 4.17 %
- 1 bp Deletion Probability: 13.32 %
- Deletion Frequency: 88.11 %
- 1 bp Frameshift Frequency: 36.35 %
- 2 bp Frameshift Frequency: 40.02 %
- Frameshift Frequency: 76.01 %

## Contact
If you encounter any issues or would like to give any feedback, feel free to leave a [GitHub issue](https://github.com/vli31/CROTON/issues).

[Back to Top](#sec0)