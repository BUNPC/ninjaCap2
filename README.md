# ninjaCap2 Repo: BRIEF DESCRIPTION
This repository contains the complete framework necessary to generate ninjaCaps. 
This help file explains how to set up and run the ninjacap code using matlab and the blender-python integration. 
Last major version: Feb 2021.

Inputs to this framework are
1. a "probe.SD" file that describes the geometrical setup (as documented in the Homer3 / AtlasViewer toolbox)
1. the desired headcircumference in cm.

Outputs of this code are four .stl files for 3D printing with ninjaflex.

# HOW TO SET EVERYTHING UP

## 1. Install Matlab and the AtlasViewer (Homer3) toolbox for Matlab
The latest Matlab version that the ninjaCap code was successfully tested with is R2019b.
Install the latest version of AtlasViewer and Homer3 (sourcecode, not the executable) for Matlab. You can get them here:	https://github.com/BUNPC/AtlasViewer
https://github.com/BUNPC/Homer3
Don't forget to install the toolbox by running "setpaths" (see AtlasViewer/Homer3 documentation).

## 2. Get the latest version of ninjaCap2 code
you are already here: make a local copy of this repository (https://github.com/neuluce/ninjaCap2).

## 3. Install Blender 2.82 64Bit 
Please note that the ninjacap2 python script may or may not be compatible with older/other versions, 
as blender python commands might be subject to changes over time.
You can get it here:		https://www.blender.org/download/releases/2-82/

## 4. Add necessary Python modules to Blender's Python
Blender uses its own Python. Blender 2.80 upwards uses Python 3.7.
For this reason, you need to install necessary modules (scipy) explicitly into Blender's Python.
If you have your own Python distribution + modules such as Anaconda installed,
blender's Python will not natively find those modules. 
If you manage modules with pip, make sure they are installed into Blender's Python!
In the instructions below, < > defines commands to run in the console and
paths are example paths that need to be adapted to the directories in your installation.
1. Start your Python CMD shell (e.g. Anaconda Prompt) or the Windows CMD console as *administrator*
1. move to the blender\python folder, (e.g. C:\Program Files\Blender Foundation\Blender\2.82\python\bin)
1. from this folder, run *<python -m ensurepip>*. This will link and install pip for python module management. Please note that Blender 2.81 and following versions might not need you to run the ensurepip command. If you get an error when running this command in version(s) 2.81+ just skip this step.
1. update pip and pip wheel by running *<python -m pip install --upgrade pip>* and
	*<python -m pip install --upgrade pip setuptools wheel>*
1. install the **scipy module** using pip. specify the target folder (--target), otherwise the modules might end up in your
	OS python folders and will not be found by Blender. E.g. *<	python -m pip install scipy --target="C:\Program Files\Blender Foundation\Blender\2.80\python\lib">* 
  (adapt directory to your folder structure accordingly).

# PREPARING FILES FOR A CAP BUILD

## 1. Create your "probe.SD" using the newest AtlasViewer GUI (ver. May2020 or newer). 
Follow the typical process to generate a probe with the AtlasViewer SDgui. Then, for each source, detector AND dummy optode, select the Grommet Type of that optode. The Grommet Type is an identifier that is used in the cap generation to place the desired element (holder, grommet, …) at the position of the corresponding optode. A brief instruction and list of currently supported grommets/holders, including their IDs, is provided [here](https://github.com/neuluce/ninjaCap/blob/master/docu/grommet_lookup.pdf). Selecting the type “#NONE” will skip the placement of a holder/ an element for this optode. You can assign any Grommet Type identifier to any optode (source, detector, and dummy optodes)!

![](https://github.com/neuluce/ninjaCap/blob/master/docu/grommetType.png)

Please note for Short Separation (SS) channels: Depending on the type of grommet you choose, short separation channels (holes for a SS fiber) are typically already included. In your SD probe layout, you might have added dedicated SS optodes. If you did, just give those a “#NONE” ID to skip the placement of a dedicated element at the SS optode location.

Please make sure you follow the keep out guidelines shown [here](https://github.com/neuluce/ninjaCap/blob/master/docu/guidelines_keepOutRegions): If you place grommets in the keep out region, they will overlap with the ear slits.
![](https://github.com/neuluce/ninjaCap/blob/master/docu/keepout.png)

# HOW TO BUILD A CAP 

## 1. Provide all relevant files:
1. place the "probe.SD" file in a new directory.
1. the following **optional** steps are available, but are not necessary by default:
    1. (optional) provide your own 'atlasViewer.mat' headmodel in the root '.../ninjaCap/' directory
    1. (optional) provide your custom side outline .svg file in the '.../ninjaCap/svg/' directory. -> change the filename in 'ninjaCap/matfun/cap.m' accordingly, if you did not use the default names 'side_piece1.svg' and 'side_piece2.svg'.
    1. (optional, advanced) provide your custom .STL elements to be placed with corresponding grommet IDs. These elements can be grommets, straps, cases, anything. Place them in '.../ninjaCap/stl/elements/<identifier>/'. For this to work, you need to follow the ID-/file-/folder- naming conventions outlined [here](https://github.com/neuluce/ninjaCap/blob/master/docu/grommet_lookup.pdf).
	
## 2. Generate the ninjaCap 
1. Open matlab, go to AtlasViewer directory and add AtlasViewer directory to matlab path by running the command 'addpath(genpath('.'))' in the matlab command line. 

2. Go to ninjaCap2 directory and add ninjaCap2 directory to matlab path by running the command 'addpath(genpath('.'))' in the matlab command line.

3. Now go to the "probe.SD" directory and run the function '* generateNinjaCap.m*'. The only argument for this function is the desired head circumference in cm. To generate a cap with 56cm headcircumference, run ' generateNinjaCap(56)'. If you do not provide an argument, 56cm is assumed per default. Watch your cap being built. By the end of the process, matlab open a Blender3D workspace.

![](https://github.com/neuluce/ninjaCap/blob/master/docu/matlab.PNG)
4. In the Blender workspace, the python script 'nynja_blender_exec.py' (ninjaCap root directory) should already be loaded per default on the right hand side. If for any case it isnt or you have your own Blender workspace, open the script window and load the script. Run the script by pressing ALT+P or clicking on 'Text -> Run Script' in the script window.

![](https://github.com/neuluce/ninjaCap/blob/master/docu/blender.png)

When blender has done its work, close it. This will also terminate the running matlab script.

## 3. Get and print the ninjaCap files
The fully generated and assembled ninjaCap panels are saved in the '.../ninjaCap2/print/' subdirectory.
You should find four .STL files there: 'sideLeft_proc.STL', 'sideRight_proc.STL', 'top1_proc.STL' and 'top2_proc_STL'.
These files can now be 3D printed. For this, use ninjaflex material with the CURA LulzBot Taz 6 Aerostruder and the printing profile 'ninjaCap_profile.curaprofile' that is also provided in the '.../ninjaCap2/print/' subdirectory.
If you want to view the files before moving on to the 3D printing step, you can open them with any STL viewer, such as for example [Meshlab](https://sourceforge.net/projects/meshlab/).

![](https://github.com/neuluce/ninjaCap/blob/master/docu/cap.PNG)


## 4. Finish the cap
With all the panels printed, you now need to assemble the cap, e.g. with an ultrasonic welder. Done.

# DIRECTORY AND FILE STRUCTURE
The ninjaCap2 repository has all main scripts and workspaces that need to be executed in the root directory. The subdirectories provide documentation, auxiliary functions, templates, user inputs, generated files, and so forth. This structure is documented in the following:

![](https://github.com/neuluce/ninjaCap/blob/master/docu/dirstruct.PNG)

* /docu/ - contains images and grommet ID manual and documentation
* /extfun/ - contains external matlab scripts, e.g. from matlab central, including the corresponding license files
* /fw/ - a default AtlasViewer directory (**ignore**)
* /imagerecon/ - a default AtlasViewer directory (**ignore**)
* /scripts/ - contains all ninjaCap matlab scripts and functions called by 'generate.m' and its subroutines.
* /print/ - contains the results of the ninjaCap generation process. You find your .STL files for 3D printing, and the 3D printing profile in here.
* /stl/ - contains (subfolders with) the cap elements that are generated by the matlab routines, as well as element templates (grommets, chinstraps, etc..). The Blender Python script assembles the cap using these elements.
* /svg/ - contains .svg files for the outline of the sidepieces. 
* /template/ - contains templates for the projection of 3D AtlasViewer probe coordinates relative to the EEG 10-10 system onto 2D cap panels.
* / - contains the necessary files for cap generation that the user needs to provide: 'probe.SD'.



