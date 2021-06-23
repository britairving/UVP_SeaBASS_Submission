# UVP_submission_formatting
https://github.com/britairving/UVP_submission_formatting/

This repository provides scripts to reformat UVP data downloaded from Ecotaxa for submission to data archival websites SeaBASS and BCO-DMO.

The Underwater Video Profiler (UVP) is designed for the quantification of particles and of large zooplankton in the water column. The UVP5 software acquires and processes images in real time. During the downcast, The UVP5 captures images of a known volume of water.  Particles are detected and analyzed to get size and grey level for each region of interest. Image post processing and metadata creation is accomplished with the Zooprocess software. Images and their associated metadata have been uploaded to the Ecotaxa website (http://ecotaxa.obs-vlfr.fr/) which serves as a tool for taxonomic annotation of zooplankton and classification of non-living particles, using machine learning and human verification, as well as a repository for all globally collected UVP data.

***
# Setup
See Setup wiki page for details on how to organize your data and configuration scripts to utilize this toolbox, and how to export data in the appropriate format from Ecotaxa's Image and Particle Modules.

https://github.com/britairving/UVP_submission_formatting/wiki/Setup
***
# Table of contents
#### Write_BCO_DMO_UVP_par.m
* User defines _cruiseid_, project directory (_projectdir_), and UVP PAR file (_odv_rfile_)
* Reads in metadata for the project from fullfile(_projectdir_, [_cruiseid_ '_UVP_metadata.m'])
* Reads in detailed ODV PAR file exported in ODV format from [Ecotaxa's Particle Module](https://ecotaxa.obs-vlfr.fr/part/)
* Format fields for BCO-DMO (since BCO-DMO does not have specific formatting requirements, this just entails renaming variables for clarity)
* Writes parameter descriptions to excel file _UVP_PAR_ParameterDescriptions.xlsx_
* Format header text (_hdr_text_) based on project metadata from {_cruiseid_}__UVP_metadata.m_ 
* Write formatted data with header to a file in _projectdir_

#### Write_BCO_DMO_UVP_zoo.m
* User defines _cruiseid_, project directory (_projectdir_), and UVP ZOO file (_odv_rfile_)
* Optionally, the user can define _validated_profiles_, a text file listing fully validated profiles. If defined, only those profiles will be written to the file.
* Reads in metadata for the project from fullfile(_projectdir_, [_cruiseid_ '_UVP_metadata.m'])
* Reads in detailed ODV ZOO file exported in ODV format from [Ecotaxa's Particle Module](https://ecotaxa.obs-vlfr.fr/part/)
* Format fields for BCO-DMO and remove unnecessary columns
* Query [WoRMs](https://www.marinespecies.org/) to match taxonomic names to scientificNames and scientificNameIDs 
* Writes parameter descriptions to excel file _UVP_ZOO_ParameterDescriptions.xlsx_
* Format header text (_hdr_text_) based on project metadata from {_cruiseid_}__UVP_metadata.m_ 
* Write formatted data with header to a file in _projectdir_

#### Write_SEABASS_Level1b_UVP_zoo.m
Level 1b definition: Individual level counts with automatic (including interpretation of class scores or probabilities ) and manual classifications, and biovolume and size parameters for each region of interest (ROI). A ROI is defined as a rectangular subset of pixels in a given image. The submission of a Level 1b data table for plankton and other particle observations to SeaBASS must include morphological information for each ROI and must be accompanied by documents that include relevant metadata and processing information.

* Read individual particle and plankton level manual classifications exported from the Ecotaxa Image module in D.O.I. export (tsv) format
* Read in metadata from [cruiseid]_UVP_metadata.m
* Read sb fields from SeaBASS_define_taxonomic_Level1b_fields.m
* Query WoRMS database for taxa scientificNames and scientificNameIDs
* Write YAML-formatted namespace file defining non-conforming ROI's
* Calculate sb fields from definitions
* Build SeaBASS formatted metadata header 
* write level1 sb file
* Write Assessed_id_list_[cruiseid].csv file

#### Write_SEABASS_Level2_UVP_par.m
* User defines _cruiseid_, project directory (_projectdir_), and UVP PAR file to read in (_odv_rfile_), SeaBASS filename to write to (_sb_wfile_)
* Optionally, the user can define _r2r_elog_, an excel sheet that lists [R2R](https://www.rvdata.us/) events with casts
* Reads in metadata for the project from fullfile(_projectdir_, [_cruiseid_ '_UVP_metadata.m'])
* Reads in detailed ODV PAR file exported in ODV format from [Ecotaxa's Particle Module](https://ecotaxa.obs-vlfr.fr/part/)
* Format fields and units for [SeaBASS](https://seabass.gsfc.nasa.gov/wiki/Data_Submission) and remove unnecessary columns
* Write formatted SeaBASS header based on project metadata from {_cruiseid_}__UVP_metadata.m_ 
* Writes a differential file, with fields PSD_DNSD, PSD_DVSD, and a non-differential file with files PSD_NSD and PSD_VSD. 

#### Write_SEABASS_Level2_UVP_zoo.m
Level 2 definition: abundance and biovolumes
* User defines _cruiseid_, project directory (_projectdir_), and UVP ZOO file to read in (_odv_rfile_), SeaBASS filename to write to (_sb_wfile_)
* Optionally, the user can define _r2r_elog_, an excel sheet that lists [R2R](https://www.rvdata.us/) events with casts
* Optionally, set _cfg.limit_to_taxa_ to a cell array of taxonomic names you're interested in
* Reads in metadata for the project from fullfile(_projectdir_, [_cruiseid_ '_UVP_metadata.m'])
* Reads in detailed ODV ZOO file exported in ODV format from [Ecotaxa's Particle Module](https://ecotaxa.obs-vlfr.fr/part/)
* Format fields and units for [SeaBASS](https://seabass.gsfc.nasa.gov/wiki/Data_Submission) and remove unnecessary columns
* Query [WoRMs](https://www.marinespecies.org/) to match taxonomic names to scientificNames and scientificNameIDs 
* Write formatted SeaBASS header based on project metadata from {_cruiseid_}__UVP_metadata.m_ 

#### WoRMS_AphiaID_taxa_match.m
* Calls the WoRMS webservice from Matlab to match taxonomic name with WoRMS AphiaID
* Utilizes [WoRMS webservice from Matlab](https://www.marinespecies.org/aphia.php?p=webservice&type=matlab)

***
# Repositories 

SeaBASS <https://seabass.gsfc.nasa.gov/#> is the SeaWiFS Bio-optical Archive and Storage System, the publicly shared archive of in situ oceanographic and atmospheric data maintained by the NASA Ocean Biology Processing Group (OBPG). 
* Example dataset: 2018 EXPORTS UVP Data <https://seabass.gsfc.nasa.gov/archive/UAF/amcdonnell/EXPORTS/EXPORTSNP/archive>

BCO-DMO <https://www.bco-dmo.org/about-us> is the Biological & Chemical Oceanography Data Management Office.
* Example dataset: Taxonomic data from P16N Repeat Hydrography cruise in 2015 <https://www.bco-dmo.org/dataset/787966>
* Example dataset: Particulate data from P16N Repeat Hydrography cruise in 2015 <https://www.bco-dmo.org/dataset/787432>
***
# Citation/acknowledgement
If you use this toolbox, please cite as follows:
_Brita Irving, (2021), UVP_submission_formatting, GitHub repository, https://github.com/britairving/UVP_submission_formatting/_
***
