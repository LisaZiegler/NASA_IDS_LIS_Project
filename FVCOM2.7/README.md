## About this Repository
This repository is structured in such a way that anyone who wants to use this tool can easily setup, execute and run FVCOM-ICM models in any coastal system. Each folder will contain a README.md, which will provide more detail. 

Main working folders
| Name | Description |
|------|-------------|


## **Steps to Running FVCOM on a Unix Server**

**...connecting to cbeps3**

ssh username@10.1.14.19
enter password

Setup directory environment from cbeps servers

Libraries required:
1. netcdf-3.6
## Create forcing files

Create your **makefile**

## Finally run the model

Before excuting the model, first locate the **makefile** in the FVCOM-2_7-model folder. Run that file.  

1. Create a screen so that you can run the model in the background and still be able to carry on with other things in the terminal.  
    - lziegler$ **screen -r**
    - to exit or detach screen **ctrl a+d**
2. If you dont want to set up a screen, another way to check if your model is still running is:  
    - lziegler$ **tail -f myrun.log** or **top**
    - to exit press **ctrl c**

The following command is used to excute the model

$ mpiexec -n 36 ../FVCOM-2_7-model/chesroms_HFX <file name>


