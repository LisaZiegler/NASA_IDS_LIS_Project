Implementing FVCOM in a Coastal Wetland-Estuarine Interface
-----------------------------------------------------------

The physical transport of water, salinity, and temperature in ICM are driven by the hydrodynamic model, FVCOM (Chen et al., 2003; Kim & Khangaonkar, 2012). FVCOM utilizes a terrain following unstructured grid, making it ideally suited for complex coastal bathymetry. In addition, the wetting and drying treatment within FVCOM makes it well suited for applications in the intertidal range (Chen et al., 2008). ICM was coupled to FVCOM for research in the Salish Sea (Khangaonkar et al., 2017; Kim & Khangaonkar, 2012). We have adapted the same modeling system used for work in Chesapeake Bay, MD, USA (Clark et al., 2020) and applied it to Long Island Sound estuary and its major freshwater tributaries. Our scale implementation of FVCOM utilizes 10 sigma layers in the vertical domain, with 14,169 elements and 7,430 nodes in a meso-haline segment of the second major river channel leading into the estuary, the Housatonic River (Figure 1a). The model resolution increases from ~700 m in the main stem of Long Island Sound to ~30 m in the marsh region

[put code to display figure here...]

The main goal of this project is to implement this FVCOM model across the entire Long Island Sound estuary, inclduing its tidal marshes. 

ICM is offline coupled to FVCOM in a system called the Unstructured Biogeochemical Model (Kim & Khangaonkar, 2012). The physical variables relating to water advection, diffusion, temperature, and salinity are first calculated independently of the biogeochemistry. Once a satisfactory hydrodynamic model solution is attained, the solution is stored and then used to drive the biogeochemical kinetic formulations within ICM. Details of all the formulations related to reaction kinetics for each of the ICM water column biogeochemical constituents (following Cerco & Noel, 2017), including organic carbon, organic nitrogen, inorganic nitrogen, underwater (UV) ultraviolet and visible light, phytoplankton, and dissolved oxygen are referenced in the Supplementary information by Clark et al., 2020 [insert citation to paper and supplementary information]. Model parameters are based on the Rhode River model and adapted where necessary for the Housatonic River system (refer to Clark et al., 2020 for detailed description of model paramters and setup).

[Insert Aim] 
[Insert questions]

Outline the structure of the folder and its contents. Providing a detailed discription with the goal in mind for anyone reading this can easily setup, execute and run FVCOM-ICM models in any coastal system.

**...connecting to cbeps3**

ssh username@10.1.14.19
enter password

Main working folders
1. FVCOM-2_7-model [**sourcecode**]
2. Housatonic      [**executablefoler**-contains input folder and run.dat excecutable file]

## **Steps to Running FVCOM on a Unix Server**

Setup directory environment from cbeps servers

Libraries required:
1. netcdf-3.6
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

$ mpiexec -n 36 ../FVCOM-2_7-model/chesroms_HFX <file name>


