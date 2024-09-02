
# EEG and MEG Methods Course Repository

Welcome to the EEG and MEG Methods course repository! This guide provides installation instructions for the software tools you'll need to follow along with the course. Depending on your preferred programming environment, choose either MATLAB or Python. If you feel savvy, you can install both.

## MATLAB Users

### 1. MATLAB Installation
If you don't have MATLAB installed, follow these steps:
1. Visit the [MathWorks website](https://www.mathworks.com/downloads/).
2. Download and install MATLAB following the provided instructions.

### 2. FieldTrip Toolbox Installation
FieldTrip is a MATLAB toolbox for MEG and EEG analysis.
1. Clone the FieldTrip repository:
    ```bash
    git clone https://github.com/fieldtrip/fieldtrip.git
    ```
2. Add FieldTrip to your MATLAB path:
    ```matlab
    addpath('/path/to/fieldtrip')
    ```
    Replace `/path/to/fieldtrip` with the actual path where you cloned the repository.

3. Initialize the FieldTrip defaults:
    ```matlab
    ft_defaults
    ```

4. (Optional) Add FieldTrip to your MATLAB startup file so it's automatically loaded each session:
    ```matlab
    edit startup.m
    ```
    Add the following line to your `startup.m` file:
    ```matlab
    addpath('/path/to/fieldtrip')
    ft_defaults
    ```
    Best practice is to add fieldtrip to the start of the scripts whenever you run them.

## Python Users

### 1. Install Visual Studio Code (VSCode)
VSCode is a powerful code editor that you'll use to write and debug Python code.
1. Download VSCode from [here](https://code.visualstudio.com/).
2. Install VSCode following the provided instructions.
3. Install the Python extension for VSCode to enhance your Python coding experience:
    - Open VSCode.
    - Go to the Extensions view by clicking on the Extensions icon in the Activity Bar on the side of the window.
    - Search for "Python" and install the Python extension by Microsoft.

### 2. Install Miniconda
Miniconda is a minimal installer for conda, a package manager that simplifies the installation of Python packages.
1. Download Miniconda from [here](https://docs.conda.io/en/latest/miniconda.html).
2. Follow the installation instructions for your operating system.

### 3. Create a Conda Environment and Install MNE
MNE is a Python package for processing and analyzing MEG and EEG data.
1. Open your terminal or command prompt and initialize conda.
    ```bash
    conda update -n base conda -y
    conda install pip -y
    conda install ipykernel -y
    ```
2. Create a new conda environment for this course:
    ```bash
    conda create -n eegmeg python=3.9
    ```
3. Activate the environment:
    ```bash
    conda activate eegmeg
    ```
4. Install MNE and other necessary packages:
    ```bash
    conda install matplotlib seaborn pandas scikit-learn numpy
    conda install conda-forge::nibabel nilearn
    pip install mne
    pip install PyQt5
    ```

### 4. (Optional) Install Jupyter Notebook, not needed if you will use VSCode.
Jupyter Notebook is a web-based interface for running Python code in an interactive environment.
1. While your `eegmeg` environment is activated, install Jupyter Notebook:
    ```bash
    conda install jupyter
    ```
2. To launch Jupyter Notebook, simply run:
    ```bash
    jupyter notebook
    ```

## Getting Started

Once you've installed the necessary tools, you can start exploring the course materials. Make sure to activate your MATLAB path or Python environment before running any code.

If you run into any issues, feel free to check the documentation of the respective tools or reach out to the course instructor or TA for help.

Happy learning!
