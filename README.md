Running FVCOM in a Coastal Wetland-Estuarine Interface
-----------------------------------------------------

## **Steps to Running FVCOM on a Unix Server**

Setup directory environment

To successfully execute your model you first need to setup your directory. I am using the terminal to do this.

1. First make sure that you have successfully installed netcdf and libraries. 
    - brew install netcdf
2. Locate the boths to those folders and update the pathways in your .bash_profile. The reason for this is so that everytime you open your terminal it first runs through those pathways set in your **.bash_profile** environment.
    - vi .bash_profile
    - insert path to those folders for e.g. lziegler$ PATH=$PATH
    - save and exit vi by typing **:wq**
  
# Create forcing files

Create your **makefile**

# Finally run the model

Before excuting the model, first locate the **makefile** in the FVCOM-2_7-model folder. Run that file.  

1. Create a screen so that you can run the model in the background and still be able to carry on with other things in the terminal.  
    - lziegler$ **screen -r**
    - to exit or detach screen **ctrl a+d**
2. If you dont want to set up a screen, another way to check if your model is still running is:  
    - lziegler$ **tail -f myrun.log** or **top**
    - to exit press **ctrl c**

The following command is used to excute the model

$ mpiexec 36 ../FVCOM-2_7-model/chesroms_HFX <file name>


