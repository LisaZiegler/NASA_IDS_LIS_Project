## About the Input Folder

Once you have successfully set up your main working directory, in this next part the user will learn about the forcings required to run the model and gain the skills to create it these files. Examples of each forcing file can be found in the zipped filled named ***Input.zip***

## Creating the Forcing Files 

*The following files are required:*
> download and unzip
1. Toolbox
2. Input.zip
3. Data folder (all data can be found in the folder ***LIS_mainresearch*** located on Storage2 external hard drive on the IMac)

***Data products used to force the model with:***
|Variable       | Product         |
|---------------|-----------------|
|Tides          | ***TPXOv9***    |
|River Discharge| ***USGS***      |
|Temperature    | ***USGS***      | 
|Weather        | ***NARR 3-day***|
|Temp/Salt_OBC  | ***CT_DEEP***   |

> the name you used to call your ***run.dat*** file is what you will used at the beginning of each named forcing file
...that is just how FVCOM identifies and calls each forcing files

*Example*
cor.dat --> **tonic**_cor.dat

### Meshgrid information
> how to process the spatial grid information...

Before you begin creating your forcing files, the meshgrid created in *SMS* must first be converted into a format compatible for MATLAB to work with moving forward. The code is found in ***Toolbox/fvcom2_7_toolbox/make_dat_files.m***

*...note: this is the meshgrid file you will be using for the rest of your pre and post processing*

> create the following files...

|Name         |Description                              | matlab script/app used     |
|-------------|-----------------------------------------|----------------------------|
|mesh.2dm     | generated meshgrid from *SMS*           | *SMS* application on vmware|
|mesh.mat     | converted 2dm file to mat file          |*make_dat_files.m*          |
|tonic_cor.dat| contain lat and lon coordinates of nodes|*make_dat_files.m*          |
|tonic_grd.dat|
|tonic_dep.dat|

### River Boundary
> create the following files...

|Name              |Description| matlab script used|
|------------------|-----------|-------------------|
|tonic_kelp_map.dat|
|tonic_riv.dat     |

### Open Boundary
> create the following files...

|Name             |Description| matlab script used|
|-----------------|-----------|-------------------|
|tonic_obc.dat    |           | *make_dat_files.m*|
|tonic_tsobc.dat  |
|tonic_elj_obc.dat|
