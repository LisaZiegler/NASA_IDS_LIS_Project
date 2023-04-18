## About the Input Folder

Once you have successfully set up your main working directory, in this next part the user will learn about the forcings required to run the model and gain the skills to create it these files. Examples of each forcing file can be found in the zipped filled named ***Input.zip***

## Creating the Forcing Files 

*The following files are required:*
> download and unzip...
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

> the name used to call the ***run.dat*** file is what the user will insert at the beginning of each named forcing file
...that is just how FVCOM identifies and calls each forcing files

*Example:*
cor.dat --> **tonic**_cor.dat

> how to process the spatial grid information...

Before creating the forcing files, the meshgrid created in *SMS* must first be converted into a format compatible for MATLAB to work with moving forward. The code is found in ***Toolbox/fvcom2_7_toolbox/make_dat_files.m***

*...note: this meshgrid file will be using for the rest of the pre and post processing analysis*

> create the following files... make reference to the files in the ***Input.zip*** folder to gain a better idea of what the file structure looks like. Most of the files will just be copy paste with few changes made regarding grid specifications in the header of the files. 

|Name         |Description                                          | matlab script/app used      |
|-------------|-----------------------------------------------------|-----------------------------|
|mesh.2dm     | generated meshgrid from *SMS*                       | *SMS* application on vmware |
|mesh.mat     | converted 2dm file to mat file                      | *convert_2dm2mat.m*         |
|blackwater_flux_map.dat| describes the transects                   | use as is;                  |
|tonic_cor.dat| x, y, lat coordinates of grid nodes                              | created from *tonic_grd.dat*; create in terminal |
|*tonic_grd.dat| specifies the geometry of the grid elements         | *make_dat_files.m*          |
|*tonic_dep.dat| specifies the geometry x,y,z (depth) @ each element | *make_dat_files.m*          |
|tonic_el_ini.dat| obc initial conditions: temp & salt            | use as is; make changes to the header|
|*tonic_kelp_map.dat| specifies marsh grass and drag by vegetation   | *WetlandForcing.m*         |
|tonic_riv.dat     | temperature, salinity and discharge @ river boundary |*make_fvcom2_7_rivers.m* |
|tonic_bfw.dat     | bottom freshwater (groundwater)                | use as is; make changes to the header|
|tonic_obc.dat     | initial temperature and salinity @ open boundary condition defined; # of nodes, node string # | *make_dat_files.m* |
|*tonic_elj_obc.dat | tidal elevation at a given time point (hourly) at each node of the open boundary| *write_elj_obc.m* |
|tonic_its.dat     | initial conditions: temperature & salinity                 | use as is; make changes to the header|
|tonic_tsobc.dat   | 3 dimensional matrix of temperature and salinit @ open boundary       | *Write_FVCOM_tsobc.m* |
|tonic_uv_ini.dat  | wind initial conditions                            | use as is; make changes to the header|
|tonic_spg.dat     | sponge layer, bottom drag- dampen velocity near the edge    | use as is; make changes to the header|
|tonic_mc.dat      | surface wind stress--> wind_non-uniform_71875.dat | *make_weather_fvcom2_7.m*|
|tonic_mc_air.dat  | weather and heat flux variables--> HFX_non_uniform_71875.dat| *make_weather_fvcom2_7.m* |

## Take note when creating the following
### - *tonic_grd.dat & tonic_dep.dat 
In the command line:
> username$ vi mesh.2dm

Then remove E3T & ND string columns
> press **esc-key** then **shift + :**
> copy&paste **% s/E3T//g**---> enter

> copy&paste **% s/ND//g**----> enter
remove the mesh info header and footer

... save as ***tonic_grd.dat***

Copy and rename *tonic_grd.dat* to *tonic_dep.dat
> username$ cp tonic_grd.dat tonic_dep.dat 

Open ***tonic_dep.dat*** 
> username$ **vi tonic_dep.dat**

Print only the last 3 columns that define depth of the elements
In the commandline enter the following
> username$ **#elements dd** (first remove the element node numbers)

> username$ **%!awk'{print $2" "$3" "$4}'**

> username$ **wq!** (save and quit)

### - *tonic_tides.dat
Files required:
- lat_lon at each grid node (text file; can obtain this when creating the *cor.dat* forcing
- time file for the year in which the model is being run for (refer to *write_elj_obc.m*). 

1. In the terminal make sure you are in the working directory where the TPXO predict tide model is --> OTPS folder
Run the following in the command line:
> username$ **make predict_tide 

> username$ ./predict_tide -ttimefile2018 <setup.inp

*In the ***setup.inp*** file, change the following:
- Point 2: insert the name of your lat_lon file
- Point 8: provide a name for the tide output file

2. Extract the M2 tides
> username$ **cat *output_name*| awk'{print $3}' >*name_M2_tides_file**

*...note: you will use this file to create the **elj_obc.dat** forcing*

### - *tonic_kelp_map.dat

Data required:
- mesh.2dm
- polygon of marsh areas in the model domain-- obtained from NWIS inventory
- download C++/Fortran libraries from Intel in order to run the *prep_kelp_map* FORTRAN model by WLong

1. First create kelp_map using ***WetlandForcing.m*** script
2. In the terminal make sure you are in the ***prep_kelp_map*** working directory 
Run the following in the command line:
> username$ ./kelp_map_gen (executable)

You will be prompted with questions to specify:
1. mesh.2dm file path and name
2. density of vegetation per m^2, diamter (m) =  200, 0.05
3. marsh drag coefficient = 0.002
4. vegetation occupies the column (sigma layers) = 1,10




