This file is a modified version of the Ecocrop Database in which the following variables were added:

- RID is the row number as appears in the Ecocrop database in the Ecocrop package
- SID is an species ID that is similar to RID but that merges same species crops. 
- FAO_Code is the "Item Code" in FAOSTAT (without leading 0s), which corresponds to the FAO Commodity List (FCL) code
	+ http://www.fao.org/economic/ess/ess-standards/commodity/en/
	+ either as a single item or grouped (e.g. beans, dry or any nes category)
	+ Some species are included in more than one FCL, so I added the columns FAO_Code2 and FAO_Code3.
- crop indicates whether it is normally cultivated for its harvestable product or not (as described in Wikipedia)
- food, fiber, pasture, forestry, and industry indicates the main use of the plant (as described in the Uses section in Wikipedia or elsewhere). 
	+ food includes, oils for cooking, spices and plants used for tea and other stimulant beverages (coffee, mate, kava), but it does not include plants used exclusively for liquors (agave) or that are chewed as stimulants (coca, qat, areca nuts) 
	+ Fodder includes grasses that may or may not be used for animal feed.
	+ Forestry includes ornamental trees and mangrove-trees
	+ industry includes medicinal plants
- plant indicates whether the species is a land plant or not (e.g. mushrooms, other fungi and algae are ommited)
- comms are comments to report inconsistencies in the data bases (such as FCL)
