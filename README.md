## World 'potential' crop species diversity

This project aims to estimate the "potential" diversity of crop species around the world.
The potential diversity is defined as the diversity that would be obtained if all suitable crop species were planted.

Thus, the metric depends on the estimation of the suitability of every crop species in each piece of cropland. 
Two methods for suitability estimation will be considered:		
 1. The FAO - Ecocrop model (as implemented in the Recocrop package in R). 		
 2. Based on Species Distribution Models (SDM)
  a. with MaxEnt considering crop presence-absence
  b. with RandomForest (or other algorithms) using crop abundance

## How to run this code:
- First run numbered scripts in main folder starting from 00
  (0000_wd file is only useful to set the working directory when running outside .Rproject)
- Then Ecocrop scripts and SDM scripts can be run independtly, following the numbers 


### Data sources:
https://www.worldclim.org/data/index.html   
http://www.earthstat.org/      
https://www.isric.org/explore/soilgrids/soilgrids-access   
http://www.fao.org/aquastat/en/geospatial-information/global-maps-irrigated-areas/latest-version    

### Crop Abundance Data Sources:
Monfreda link
SPAM link



*Some things to be considered in future steps*		
 * Model the effect of crop demand on crop diversity: to see how much of the "diversity gap" is explained by a greater demand for a few crops.		
 	- This could be coupled with an analysis of how changes in demand would change diversity		
 	- We can also consider trade, as it has been shown that, within countries, demand diversity is greater than production diversity... but most countries can't produce everything they consume. 

	

