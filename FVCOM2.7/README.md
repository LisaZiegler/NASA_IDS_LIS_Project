## About this Repository

Main working folders
| Name | Description |
|------|-------------|

## Model domain
Our scale implementation of FVCOM utilizes 10 sigma layers in the vertical domain, with 14,169 elements and 7,430 nodes in a meso-haline segment of the second major river channel leading into the estuary, the Housatonic River. The model resolution increases from ~700 m in the main stem of Long Island Sound to ~30 m in the marsh region.

![](./github-figures/mesh_hr.jpeg)

The main goal of this project is to implement this FVCOM model across the entire Long Island Sound estuary, including its tidal marshes. Specifically focusing on assessing changes in estuarine water quality, organic matter cycling. We are currently in phase 1 of seting up the LongIsland Sound model (figure 3).

![](./github-figures/lis_grid.001.jpeg)

## **How to- Running FVCOM on a Unix Server**

**...connecting to cbeps3**

ssh username@10.1.14.19
enter password

Setup directory environment from cbeps servers

Libraries required:
1. netcdf-3.6

### Create forcing files

Create your **makefile**

### Finally run the model

..Yay you are ready to run!


1. Create a screen so that you can run the model in the background and still be able to carry on with other things in the terminal.  
    - lziegler$ **screen -r**
    - to exit or detach screen **ctrl a+d**
2. If you dont want to set up a screen, another way to check if your model is still running is:  
    - lziegler$ **tail -f myrun.log** or **top**
    - to exit press **ctrl c**

The following command is used to excute the model

$ mpiexec -n 36 ../FVCOM-2_7-model/chesroms_HFX <file name>


