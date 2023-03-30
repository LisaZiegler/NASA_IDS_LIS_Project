## About the Input Folder

Once you have successful set up your main working directory. In this next part the user will learn about the forcings required to run the model and gain the skills to create it these files. Examples of each forcing file can be found in the zipped filled named ***Input.zip***

## Creating the Forcing Files 

*The following files are required:*
1. Download the toolbox
2. Input.zip
3. Data folder (all data can be found in the folder ***LIS_mainresearch*** located on Storage2 external hard drive on the IMac)

Before you begin creating your forcing files, the meshgrid created in *SMS* must first be converted into a format compatible for MATLAB to work with moving forward. The code is found in ***Toolbox/fvcom2_7_toolbox/make_dat_files.m***

*...note: this is the meshgrid file you will be using for the rest of your pre and post processing*

These are the forcings related to the meshgrid information:

|Name         |Description                              | matlab script/app used     |
|-------------|-----------------------------------------|----------------------------|
|mesh.2dm     | generated meshgrid from *SMS*           | *SMS* application on vmware|
|mesh.mat     | converted 2dm file to mat file          |*make_dat_files.m*          |
|tonic_cor.dat| contain lat and lon coordinates of nodes|*make_dat_files.m*          |


### Physical Forcings
### River Boundary
The following files will need to be created:

|Name|Description|
|----|-----------|


### Open Boundary
The following files will need to be created:

|Name         |Description| matlab script|
|-------------|-----------|--------------|
|tonic_obc.dat|           | *make_dat_files.m*|

