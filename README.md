# Pytheas installation and basic usage notes

## General notes
Pytheas has been fully tested on Windows 10 and Linux (Ubuntu 20.04+) systems, and not fully tested but working on macOS systems. In the current release, there are two supported versions, one with a graphical interface (GUI) and one exclusively command line based (CL). The GUI version is the preferred one, but in case the user experiences issues while using it, the CL version is a good alternative albeit less user-friendly. For the purposes of inputs and outputs the scripts of the two versions are identical. Both versions require to call the component of the Pytheas workflow via command line, but the GUI version offers to the user a graphical interface to select the parameters and run the analysis. A third option is to run Pytheas through Docker, using the provided Dockerfile to build the image. 
Several instances of files and libraries needed for the execution of other scripts (either already present or created while running Pytheas) are packaged together in the same directory, please be careful if moving the scripts from their original position as it may break some dependencies. It is therefore suggested to run the scripts in their provided directories and transfer around only the output files generated.

## Pytheas Installation
### Download GitHub repository
From the [GitHub repository](https://github.com/ldascenzo/pytheas) Pytheas can be downloaded either selecting the green “Code” button and choosing “Download ZIP” or cloning locally the git repository. The first option downloads a compressed directory with all the files, including the two GUI version and CL version main working directories. To clone the git repository, install git ([Windows instructions](https://git-scm.com/download/win)), then open a shell prompt (in Windows under *Start > Command Prompt*) and navigate into the desired directory for cloning Pytheas repository. Then, within this directory run `git clone https://github.com/ldascenzo/pytheas.git`
This will create a *pytheas* directory at the chosen location, with all the files needed to run Pytheas.

### Windows 10
Download the latest [Anaconda version](https://docs.anaconda.com/anaconda/install/) based on the Windows version in usage (by April 2021 it is *Python 3.8 Windows 64-bit Graphical Installer*).
Install Anaconda leaving the standard options selected. 

Open an Anaconda Powershell Prompt (*Start menu -> Anaconda Powershell Prompt*).
Within the prompt install some Python packages needed for Pytheas execution.
```
conda install -y -c conda-forge gooey
conda install -y -c conda-forge biopython
conda install -y -c bioconda pyteomics
```

### Linux (Ubuntu 20.04)
Install the latest Python 3 version
```
$ sudo apt-get update
$ sudo apt-get install python3
```

Install pip, Python package manager to easily download, update and maintain all the libraries
```
$ sudo apt install -y python3-pip
```

Install some additional libraries that will be needed later, confirming with *Y* when needed during the installation. 
```
$ sudo apt-get install git curl libsdl2-mixer-2.0-0 libsdl2-image-2.0-0 libsdl2-2.0-0
```

Install the wxPython package, which will be used for the GUI version. 
```
$ wget  -P ~ "https://extras.wxpython.org/wxPython4/extras/linux/gtk3/ubuntu-20.04/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl"
$ pip3 install ~/wxPython-4.1.0-cp38-cp38-linux_x86_64.whl
```

Navigate to the local pytheas directory, then install all the requirements needed from the file *requirements.txt*. In case of problems with any of those, please refer to the respective documentations. 
```
$ pip3 install -r requirements.txt
```

### macOS support
The support on macOS systems is not fully tested, but the Anaconda installation described for the Windows version is the recommended way of installing Pytheas dependencies. On top of the Windows instructions, the requirements to run the GUI version are installing python.app with

```conda install -c anaconda python.app```

And calling the scripts using `pythonw` 

## Pulling Pytheas repository through Docker
### Windows 10
Install [Docker Desktop](https://hub.docker.com/editions/community/docker-ce-desktop-windows)
Open a Command Prompt (from the Start menu) and pull the Pytheas Docker repository
```
docker pull ldascenzo/pytheas:latest
```

### Linux (Ubuntu 20.04)
Follow the installation instructions for [Docker](https://docs.docker.com/engine/install/ubuntu/)
Open a Terminal prompt to pull the Pytheas Docker repository
```
docker pull ldascenzo/pytheas:latest
```

## Usage info for Linux, Windows and Docker versions
### Windows 10
Using an *Anaconda Powershell Prompt*, navigate to the local pytheas directory. If for example the directory is in C:\Users\MainUser\pytheas the command will be
`cd C:\Users\MainUser\pytheas`

Navigate inside the directory of the chosen version (by default the GUI version) and run the scripts inside the respective sub-directories by using 
```python pytheas_<script_name>.py```

Select the options and parameters using the provided graphical interface. Refer to the Fig.1 of the main text and to the manual for more details on the order of the workflow and the usage/output of the single scripts.
### Linux (Ubuntu 20.04)
Navigate within the Pytheas version directory (GUI or CL) and launch the scripts with 
```python3 pytheas_<script_name>.py```

Note that python3 instead of python may not be necessary if the default Python installation is set to be the 3. Refer to the Fig.1 of the main text and to the manual for more details on the order of the workflow and the usage/output of the single scripts.

### Docker container
Check the name of the newly pulled container under the column “NAMES”
```docker container ls```

Open a Bash shell within the container using the name you just recovered 
```docker exec -it <container name> /bin/bash```

Navigate within the main Pytheas directory, in the command line version
```# cd pytheas/CL_version```

From there run Pytheas scripts following the normal instructions for the CL version
To transfer files between the Docker container and the host, use the following
```docker cp <container name>:/file/path/within/container /host/path/target```

