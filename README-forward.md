# Analysis of GROMACS 2020 simulations with Jupyter notebooks
## Introduction
We have prepared jupyter notebooks which allow a more straightforward analysis of GROMACS simulations.
Currently the notebooks are based on running and collecting the output from GROMACS programs, but
we also use the MDAnalysis toolkit for some analyses. Since a browser is required we recommend two methods of using the notebooks to analyse data from CINECA Marconi M100:
1. Install the python software on your own PC, copy the data down from M100 and run jupyter locally.
2. Use the forward utility (link below) to run jupyter on M100 but visualizing the output on your own pc (via the browser).


In the second case you should install the python software via a virtual environment in your home directory on M100 (see below). For local installations it is not essential to use a virtual environment, but it is recommended. 

## Requirements
1. python3
2. web browser
3. notebook file and library (`EX4Cutils.py`)
4. GROMACS 2020.x or higher
5. Machine installed scipy module.

## Setup for analyses of data on M100

### Local workstation installations for visualizing notebook output
You need native Linux or Linux running under WSL (Windows).
You first need to create and SSH key and copy it to M100.
Then install the forward software so you can visualise the notebook on a PC browser.
```
ssh keygen -t rsa 
ssh-copyid username@login.m100.cineca.it     # username=M100 username
git clone  https://gitlab.hpc.cineca.it/mrorro00/forward2.git
```
Check to make sure you can log on to M100 without a password.

### Using a pre-installed environment on M100
Now logon to M100 and clone this repository, copying files to home directory
```
git clone https://gitlab.hpc.cineca.it/training/molecular-dynamics.git
cd molecular-dynamics/Tutorials/Cluster-Analysis
cp analysis-venv EX4Cutils.py  single-analysis.ipynb $HOME
```

To avoid installing the python software we have prepared a pre-installed python environment - ```analysis-venv```.
*Caution this will only work for CINECA's M100 system*.

Although pre-compiled you still need to load some modules on M100.
To check the virtual environment is correctly installed:
```
module load autoload python
module load profile/deeplrn autoload scipy
source analysis-venv/bin/activate
```

### Manual M100 Setup of virtual environment
***Not needed unless pre-compiled environment does not work***

Log onto M100, load modules and set up virtual environment-make sure you have Python3: 
```
module load autoload gromacs 
module load python
module load profile/deeplrn autoload scipy
python -m venv analysis
source my_venv/bin/activate
```
Now install python libraries as described in the Library install section below.

## Analysis of data on local workstation
Make sure you have python3 and GROMACS 2020 installed.
We recommend a virtual environment also for local installations:
```
python -m venv analysis
source analysis/bin/activate
```

### Library install
*Warning: this may take some time, but only needs to be done once.*
*skip if using the pre-installed environment*
```
pip3 install matplotlib pandas
pip3 install MDAnalysis
pip3 install numpy # unless loaded via a module
pip3 install notebook
```
If you have difficulty installing Pandas then try:
```
pip3 install ipykernel --upgrade
python3 -m ipykernel install --user
```
Difficulties with matplotlib can be resolved by  loading cython:
```
pip install cython
pip install matplotlib
```

For visualization with nglview (not currently working on M100):
```
pip3 install nglview
pip3 install ipywidgets
jupyter nbextension enable --py widgetsnbextension
```
If not using a virtualenv you may need to the last command with sudo (for local installations).

## Analysis
### Running the analysis notebook locally
First make sure you that you have the `EX4Cutils.py` library file in the same directory from where you run
jupyter. Then launch:

```
jupyter notebook
```
A browser should open automatically, otherwise copy and paste the link provided into an open one.
Within the browser locate and open the ```single-analysis.ipynb``` file.

### Running the analysis on M100 using the forward utility
The forward utility by Marco Rorro launches jupyter on the M100 but tunnels the output to your own pc so you can
visualise the analysis on a browser. It uses the SLURM batch system so you must first create the batch file. We have provided one so do the following.
1. Copy it from the course repository or from M100 onto your workstation:
```
scp username@login.m100.cineca.it:molecular-dynamics/Tutorials/Cluster-analysis/jupyter-notebook .
```
Or use ```git clone``` to copy the tutorial directory on your local workstation.

2. Modify ```jupyter-notebook``` to correctly identify your account (```--account```) and any other job-specific options.
3. copy the file into the forward m100 directory
```
cp jupyter-notebook forward2/m100
```

4. launch jupyter on m100 using the following command:
```
cd  forward2
./forward2 -u username -r m100 jupyter-notebook
```
Where ```username``` is your M100 username.
Once you jupyter has been sucessfully launched on m100, you will be provided with link which should paste in 
a browser on your pc.

### Hints
Hints for using jupyter notebooks:
- cell types can be python code or markdown - this can be changed from the menu
- <enter> inserts a carriage return, while <shift><enter> executes the cell if containing python code or renders
the cell content, if containing markdown
- cells containing % will run unix commands in the shell
- output is plotted in line with matplotlib


### Performing the analysis
1. Go to "Define files" and set the locations of the files to analyse.
For gromacs analysis you need for each simulation:
    - a .tpr file
    - a pdb file 
    - an .xtc (with pbc removed)
    - possible index files, for holo complexes. 

We recommend you use the first frame of the xtc file to obtain the pdb.
```
gmx trjconv -f file.xtc -dump 0 -o file.pdb
```

In the current version of EX4Cutils you need an index. This can be created from GROMACS as,
```gmx make_ndx -f topol.tpr -o index.ndx```
Choose 'q' to exit and save.


2. RMSD
This has been implemented as an MDAnalysis command (so you wont see GROMACS output in the terminal window).
From the RMSD output you can estimate the initial time (given by the *start* variables) to discard by the adjusting the ```xlim[]``` parameter in the 
plot.

3. PCA and Clustering
This is done by GROMACS. If output from the notebook is not as expected , then check the terminal window to see if GROMACS indicated an error. Note that the cluster and PCA analysis can be very time consuming.

### Saving the analysis
Apart from the notebook it is also possible to save the analysis in other ways:
- pdf files (easiest way is via the browser)
- html
```bash
jupyter nbconvert --execute --to html notebook.ipynb
```
4. Suggested workflow
   1. copy notebook file and EXC4utils.py to directory with trajectory
   2. Create pdb file with trjconv and options -dump 0 -o file.pdb
   3. Create index file with gmx make_ndx
   4. Run jupyter notebook and navigate directly to the directory.

 
