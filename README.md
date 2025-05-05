# sim_jtb2012
# Oocyte simulations supporting JTB 2012 Paper

This project is a graphical user interface (GUI) front end to the cellular modelling simulaions discussed in the paper "Somersalo, E., et al., A reaction–diffusion model of CO2 inﬂux into an oocyte. J. Theor. Biol. (2012), http://dx.doi.org/10.1016/j.jtbi.2012.06.016".  The goal: to make it easier for users to run simulations and visualize the results; and to allow for batch running of large numbers of simulations.  The orignal simulation code was written in Matlab.  You can find the original files in the folder OriginalMatlabCode.  See the License file for authorship, licensing information.

This project has been written in Python, with hooks to Matlab, running modified versions of some the original Matlab files, to take in parameters from the GUI and retrieve the results.  As of the time of this was written, Matlab is required to run the code.  Future versions are planned for running with Octave and an all python version is also planned.


## Installation:
    Linux:
    1. Pre-requisites:
        Matlab ( last tested with 2024b  )
        Python ( last tested with 3.12.3 )
    2. Create virtual env and install the required packages
        from the same directory where you clone this git repo
        > python3 -m venv venv
        > source venv/bin/activate
        > pip install matplotlib numpy pandas pycairo scipy -f https://extras.wxpython.org/wxPython4/extras/linux/gtk3/ubuntu-24.04/ wxPython
    3. Install matlab engine
        find path of your Matlab install with command "matlabroot" from Matlab command window, cd to that folder, then to extern/engines/python  for example if Matlab is installed in your home folder:
        > cd ~/MATLAB/R2024b/extern/engines/python/
        > pip install .

    Windows: TODO: Instructions to follow at a later date.


## Usage:
    At this point you should be able to run the program with
  > python mgui.py
  -h or --help will give you many command line options which I find to be the most efficient manner to run the program.

The GUI is designed around the idea that it should be easy for the users to modify parameters fed to the simulation code.  With future expansion in mind it will support other simulation types, soon AJP and RBC will be added in the next releases.

In the top left you can choose a simulation type, a paper and a figure, which will fill in the default starting paramters for the selected simulation.  Selecting "JTB", the paper, and the specific figure will load a "Parameter" set from a Param file that contains the defaults needed to accurately recreate the figures in the paper.  For example you can run Figure "Figure_3_4_5", this simulation is used to create Figures 3, 4 and 5 from the paper.  When the simulation finishes you can push the button for "Fig 3" "Fig 4" or "Fig 5" to generate that figure.  You can modify any paramaters you wish, save that as a custom parameter set, load a custom set and load data file that you have created previously.  Most simply, hit "Run Simulation", answer folder prompts and wait for the simulation to finish, and hit the Figure button you wish to see.

Default parameters for the figures in the paper are stored in the "Params" folder.  For custom simulations, go to Params/Sim_1__JTB/Paper_2__Nonspecific and create a new Param file containing your custom values.  They are python files, follow the same format (TODO: better document this).

Figures code is stored in teh "Figures" folder.  Notice that Params and Figures follow the same hierarchy.  (TODO: better document this)

