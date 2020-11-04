# Information about the data used to force the open and river boundaries for both the Housatonic and Long Island Sound Models

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

- Won't be representing the marsh in the LIS model. 

## Open Boundary

- For the Housatonic Model a monthly average of temperature, salinity was derived from 9 CT DEEP stations (LISICOS) in the sound,
this is because the WOA (World Ocean Atlas) data would overestimated temperature and salinity (create bias)
- For the LIS Model the WOA model was used to obtain temperature and salinity. Tides m2 (TPXOv9 model). Weather from NARR (3hr composites)



