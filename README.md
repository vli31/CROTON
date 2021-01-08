<a id='sec0'></a>
# CROTON
Please follow these instructions to host CROTON on your local computer. Reference this link for more information about CROTON and its capabilities: [Link to Paper](https://drive.google.com/file/d/15eP01jxLilJnERmMnDkwn7kYkOKwOs-I/view?usp=sharing)

### Get the latest code
Clone the Github Repository. If you have a previous version, make sure to pull the latest commits:
```
git clone https://github.com/vli31/CROTON
cd CROTON
git pull
```
If you see `Already up to date` in the terminal, it means the code you have is the latest version.

## Creating an Anaconda microenvironment
To run CROTON, it is easiest to use a conda microenvironment. If you do not have Anaconda, you can download it at this link: https://www.anaconda.com/products/individual. 

The environment used for CROTON is available in this repository in the file 'croton.yml'. It is recommended that you use the same versions of Django (3.0.3), Keras (2.2.4) and TensorFlow (1.14.0).
```
conda env create -f croton.yml
conda activate croton
```

## Setting up Django 
To set up django make sure you are in the directory 'CROTON' and have activated the microenvironment 'croton', then run:
```
python manage.py runserver
```
Now CROTON should be avilable at a local host like http://127.0.0.1:8000/

The CROTON web interface should look like this:

<img pull-left src="croton_screenshot.png">

## Contact
If you encounter any issues or would like to give any feedback, feel free to leave a [GitHub issue](https://github.com/vli31/CROTON/issues).

[Back to Top](#sec0)