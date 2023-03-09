## Implementing FVCOM in a Coastal Wetland-Estuarine System

The physical transport of water, salinity, and temperature in ICM are driven by the hydrodynamic model, FVCOM (Chen et al., 2003; Kim & Khangaonkar, 2012). FVCOM utilizes a terrain following unstructured grid, making it ideally suited for complex coastal bathymetry. In addition, the wetting and drying treatment within FVCOM makes it well suited for applications in the intertidal range (Chen et al., 2008). ICM was coupled to FVCOM for research in the Salish Sea (Khangaonkar et al., 2017; Kim & Khangaonkar, 2012). We have adapted the same modeling system used for work in Chesapeake Bay, MD, USA (Clark et al., 2020) and applied it to Long Island Sound estuary and its major freshwater tributaries. 

![](./github-figures/bathy_bw.001.jpeg)

## About FVCOM-ICM modeling tool

ICM is offline coupled to FVCOM in a system called the Unstructured Biogeochemical Model (Kim & Khangaonkar, 2012). The physical variables relating to water advection, diffusion, temperature, and salinity are first calculated independently of the biogeochemistry. Once a satisfactory hydrodynamic model solution is attained, the solution is stored and then used to drive the biogeochemical kinetic formulations within ICM. Details of all the formulations related to reaction kinetics for each of the ICM water column biogeochemical constituents (following Cerco & Noel, 2017), including organic carbon, organic nitrogen, inorganic nitrogen, underwater (UV) ultraviolet and visible light, phytoplankton, and dissolved oxygen are referenced in the Supplementary information by Clark et al., 2020 [insert citation to paper and supplementary information]. Model parameters are based on the Rhode River model and adapted where necessary for the Housatonic River system (refer to Clark et al., 2020 for detailed description of model paramters and setup).

The model allows exploration of the role tidal marshes play in estuary carbon and nitrogen cycling by simulating physical and biogeochemical processes at a high spatial and temporal resolution. By numerically removing the tidal marsh from the model, a direct estimate is made of the total DOC that can be attributed to the presence of the tidal marsh. This has implications for how much DOC in estuarine organic carbon budgets can be considered of tidal marsh origin. It should be noted that the marsh area itself does not contain a dynamic plant community and the associated biogeochemical effects. The lack of dynamic plant growth and senescence means that the TOC budget within the marsh cannot be estimated; the marsh acts more as a source or sink for organic matter rather than a dynamic part of the ecosystem.  

## About this Repository
This repository is structured in such a way that anyone who wants to use this tool can easily setup, execute and run FVCOM-ICM models in any coastal system. Each folder will contain a README.md, which will provide more detail. 

*These are the main working folders you will be using:*

| Name | Description |
|------|-------------|
| FVCOM2.7- physics| includes: sourcecode, executable folder (inputs and run.dat executable file)|
| ICM- biogeochemistry  | includes: sourcecode, executable folder (inputs and run.dat executable file)|
| Toolbox| contains matlab scripts used to make model forcing files and plot model output|
| INSTALL_modules| these modules need to be installed before running both models. The sourcecode calles on these libraries during model computation|

------
## Information about the data used to force the open and river boundaries for both the Housatonic and Long Island Sound Models

## River Boundary Forcing

- The Thames and Connecticut Rivers make up 80% of the freshwater inflow into the sound. With the Connecticut contributing 72%,
the largest freshwater input source

These are the princple rivers flowing into the LIS

1. Hudson River 
2. Housatonic River (combines Naugatuck and Housatonic Rivers)
3. Quinnipiac River
4. Connecticut River
5. Thames River (combines Quinebaug, the Shetucket, and the Yantic Rivers)

USGS gauge stations used

1. Hudson @ Green Island, Troy Dam -- 01358000 (discharge); 01359139 (daily mean temperature 2014-2016)
2. Naugatuck @ Beacon Falls -- 01208500 (discharge), Housatonic @ Stevenson -- 01205500 (discharge), Milford -- 01200600 (temperature)
3. Quinnipiac River @ Wallingford -- 01196500 (discharge and temperature)
4. Connecticut @ Thompsonville -- 01127500 (discharge and temperature)
5. Yantic River @ Yantic -- 01127500 (discharge, Quinebaug River at Jewett City -- 01127000 (discharge), Shetucket River near Willimantic -- 01122500 (discharge)

For the Housatonic Model domain a sum of the freshwater discharge was obtained by combining the Naugatuck and the Housatonic Rivers. 
This was done because there is no gauge station near the head of the Housatonic River, and thus unsure of what that flow might be. 
Connecticut will be used for the Thames, assuming that the systems share similar charactertistics since they are closer to the MAB.

In terms of the LIS Model, the same will be done. All rivers are forced at the head of the river where there is no tidal influence. Since we are unsure of the exchange with the East River on the western part of the sound, the Hudson is included as a river (FVCOM computes the exchange anyway so we can validate model with 
Gay et al. 2004). All the other smaller non-point source rivers will be delt with as a shoreline, since its contribution to the freshwater flux into the sound is very small (insignificant). 

## Open Boundary

- For the Housatonic Model a monthly average of temperature, salinity was derived from 9 CT DEEP stations (LISICOS) in the sound,
this is because the WOA (World Ocean Atlas) data would overestimated temperature and salinity (create bias)
- For the LIS Model the WOA model was used to obtain temperature and salinity. Tides m2 (TPXOv9 model). Weather from NARR (3hr composites)

