## About the Input Folder

Once you have successfully set up your main working directory, in this next part the user will learn about the forcings required to run the model and gain the skills to create it these files. Examples of each forcing file can be found in the zipped filled named ***Input_2019.zip***

## Creating the Forcing Files 

*The following files are required:*
> download and unzip...
1. Toolbox
2. Input.zip
3. Data folder (all data can be found in the folder ***LIS_mainresearch*** located on Storage2 external hard drive on the IMac)

***Data products used to force the model with:***
|Variable       | Product         |
|---------------|-----------------|
|Tides          | fvcom output    |
|Temperature & Salinity | fvcom output | 
|Weather        | fvcom output    |
|River Discharge| ***USGS***      |
|Water quality  | **USGS** (riv) ***CT_DEEP*** (obc) |

> create the following files... make reference to the files in the ***Input.zip*** folder to gain a better idea of what the file structure looks like. Most of the files will just be copy paste with few changes made regarding grid specifications in the header of the files.

|Name         |Description                                          | matlab script/app used      |
|-------------|-----------------------------------------------------|-----------------------------|
|mineralisation_newDOM.3_yr_8| organic matter  mineralisation, nutrient paramterisation file       | used as is; make changes directly to file if needed |
|seddom_control.in| sediment DOM input control file; controls reactions rates of DOM and diffusion to the overlying water column | used as is; make changes directly to file if needed |
|wcdom_apqy.npt| apparent quantum yield parameterization| used as is; make changes directly to file if needed|
|wcdom_control.npt| water column DOM input control file| used as is; make changes directly to file if needed|
|wqm_algae.3_yr_7| algal paramterisation file | used as is; make changes directly to file if needed |
|wqm_kei_photoDEG_v8_1.csv| specifies water column light absorption| used as is |
|settling.10_yr_3| particulate sinking rate parameterisation | used as is; make changes directly to file if needed |
|seddom_grid.dat| sediment DOM/marsh grid definition file           | *make_seddom_gridm*|
|tonic_initial_wq_vert_restart1.dat| specify initial conditions| *write_ICM_initial_ncdf.m*|
|tonic_met_good.dat| meteorological data (NARR) | *write_metgrid2file.m*|
|tonic_obc_wq.dat| open boundary water quality conditions timeseries over the water column| *WQM_OBC_DoCoFVCOM.m*|
|tonic_pnt_wq.dat| river point source forcing; turned off| not used|
|tonic0000_riv.dat| river water quality variables specified| *Write_DoCoFVCOM_rivers.m*|
|tonic_seddiag_2b.npt| ediment diagenesis parameter file| used as is; make changes directly to file if needed|
