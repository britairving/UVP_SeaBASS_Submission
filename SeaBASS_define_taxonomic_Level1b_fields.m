function cols = SeaBASS_define_taxonomic_Level1b_fields
%% FUNCTION SEABASS_DEFINE_TAXONOMIC_LEVEL1B_FIELDS
% Description:
%   Defines fieldnames required for Level1b SeaBASS files
%
% References:
%
%  OCB_PTWG_TM_FinalDraft_Dec2020.docx
%
%  https://seabass.gsfc.nasa.gov/wiki/plankton_and_particles
%
%  https://seabass.gsfc.nasa.gov/wiki/stdfields
%
%  Ecotaxa field definitions: https://sites.google.com/view/piqv/ecotaxa?authuser=0
%
%  level 1b definition: Individual level counts with automatic (including
%  interpretation of class scores or probabilities ) and manual
%  classifications, and biovolume and size parameters for each region of
%  interest (ROI). A ROI is defined as a rectangular subset of pixels in a
%  given image. The submission of a Level 1b data table for plankton and
%  other particle observations to SeaBASS must include morphological
%  information for each ROI and must be accompanied by documents that
%  include relevant metadata and processing information.
%
% Author:
%  Brita Irving <bkirving@alaska.edu>
%% The following is a list of field names (i.e., measurement labels) that align with Darwin Core when possible.
cols = struct();
cols.associatedMedia.description = 'A unique persistent identifier of the media associated with the occurrence. The field provides the unique imagery file name corresponding to the source of the ROI. Alternately, use this field to provide a URL pointing to a permanent landing page for the ROI image. In the latter case, instructions should be provided as comments in the header on how to construct the local file name based on the URL. If the local imagery file name cannot be constructed from the URL, then list both filenames. Use the pipe character ?|? to separate the names, and do not use spaces.';
cols.associatedMedia.requirement = 'required';
cols.associatedMedia.units       = 'none';
cols.associatedMedia.ecotaxavar  = 'img_file_name';

cols.data_provider_category_automated.description = 'A category used by the data provider to name the organism or particle for an automated classification, not necessarily a scientific name (e.g., pennate or detritus).';
cols.data_provider_category_automated.requirement = 'optional'; % (recommended but optional) 
cols.data_provider_category_automated.units       = 'none';
cols.data_provider_category_automated.ecotaxavar  = '';

cols.scientificName_automated.description = 'A scientific name from a recognized taxonomic reference database (e.g., WoRMS, AlgaeBase) at the lowest level that matches the data provider?s category for an automated classification paired to a scientificNameID. Generally, the ROI corresponds to an occurrence assigned to a single taxonomic name.'; 
cols.scientificName_automated.requirement = 'optional'; % (recommended but optional) 
cols.scientificName_automated.units       = 'none';
cols.scientificName_automated.ecotaxavar  = '';

cols.scientificNameID_automated.description = 'A life science identifier (LSID) from a recognized taxonomic reference database (e.g., WoRMS, AlgaeBase) at the lowest level that matches the data provider?s category for an automated classification. e.g., urn:lsid:marinespecies.org:taxname:233015 where ?urn:lsid? indicates the ID that is specific to life science data and is used for all files, marinespecies.org is the url for the reference database WoRMS, and the namespace ?taxname? informs the user that the following number represents a unique numerical identifier or taxon identifier in WoRMS. 233015 represents the taxon identifier (AphiaID) in WoRMS for the dinoflagellate species Karenia brevis.';
cols.scientificNameID_automated.requirement = 'required'; 
cols.scientificNameID_automated.units       = 'none';
cols.scientificNameID_automated.ecotaxavar  = '';

cols.data_provider_category_manual.description = 'A category used by the data provider to name the organism or particle for a manual identification, not necessarily a scientific name.';
cols.data_provider_category_manual.requirement = 'optional'; % (recommended but optional) 
cols.data_provider_category_manual.units       = 'none';
cols.data_provider_category_manual.ecotaxavar  = '';

cols.scientificName_manual.description = 'A scientific name from a recognized taxonomic reference database (e.g., World Register of Marine Species, AlgaeBase) at the lowest level that matches the data provider?s category, for a manual identification matched to scientificNameID. Generally, the ROI corresponds to an occurrence assigned to a single taxonomic name.';
cols.scientificName_manual.requirement = 'optional'; % (recommended but optional) 
cols.scientificName_manual.units       = 'none';
cols.scientificName_manual.ecotaxavar  = '';

cols.scientificNameID_manual.description = 'An LSID from a recognized taxonomic reference database (e.g., World Register of Marine Species, AlgaeBase) at the lowest level that matches the data provider?s category for a manual identification.';
cols.scientificNameID_manual.requirement = 'required'; %
cols.scientificNameID_manual.units       = 'none';
cols.scientificNameID_manual.ecotaxavar  = '';

cols.biovolume.description = 'Biovolume for the target detected within the ROI determined by means specified in the biovolume calculation method or protocol document.';
cols.biovolume.requirement = 'required'; %
cols.biovolume.units       = 'um^3';  
cols.biovolume.ecotaxavar  = ''; 

cols.area_cross_section.description = 'Cross-sectional area of the target detected within the ROI determined by means specified in the image processing method or protocol document.';
cols.area_cross_section.requirement = 'required'; %
cols.area_cross_section.units       = 'um^2';  
cols.area_cross_section.ecotaxavar  = '';

cols.length_representation.description = 'Representation of length of the target detected within the ROI or largest mesh size for which the target could be retained, determined by means specified in the image processing method or protocol document.';
cols.length_representation.requirement = 'required'; %
cols.length_representation.units       = 'um'; 
cols.length_representation.ecotaxavar  = '';

cols.width_representation.description = 'Representation of width of the target detected within the ROI or smallest mesh size through which the target could pass, determined by means specified in the image processing method or protocol document.';
cols.width_representation.requirement = 'required'; %
cols.width_representation.units       = 'um'; 
cols.width_representation.ecotaxavar  = '';

cols.equivalent_spherical_diameter.description = 'Equivalent spherical diameter of the target detected within the ROI determined by means specified in the image processing method or protocol document.';
cols.equivalent_spherical_diameter.requirement = 'optional'; %
cols.equivalent_spherical_diameter.units       = 'um'; 
cols.equivalent_spherical_diameter.ecotaxavar  = '';

cols.area_based_diameter.description = 'Area-based diameter of the target detected within the ROI determined by means specified in the image processing method or protocol document.';      
cols.area_based_diameter.requirement = 'optional'; %
cols.area_based_diameter.units       = 'um'; 
cols.area_based_diameter.ecotaxavar  = '';

% Fields not listed in SeaBASS wiki on plankton_and_particles (see
% comments) but important to note
cols.depth.description = 'Depth of measurement';
cols.depth.requirement = 'optional'; %
cols.depth.units       = 'm'; 
cols.depth.ecotaxavar  = 'object_depth_min';

cols.lat.description = 'Depth of measurement';
cols.lat.requirement = 'optional'; %
cols.lat.units       = 'm'; 
cols.lat.ecotaxavar  = 'object_lat';

cols.lon.description = 'Depth of measurement';
cols.lon.requirement = 'optional'; %
cols.lon.units       = 'degrees'; 
cols.lon.ecotaxavar  = 'object_lon';

cols.date.description = 'Depth of measurement';
cols.date.requirement = 'optional'; %
cols.date.units       = 'yyyymmdd'; 
cols.date.ecotaxavar  = 'object_date';

cols.time.description = 'Depth of measurement';
cols.time.requirement = 'optional'; %
cols.time.units       = 'hh:mm:ss'; 
cols.time.ecotaxavar  = 'object_time';
cols.time.convert     = true;

% Below are Required in the file header or as a field
cols.volume_sampled_ml.description = 'Original volume of sample collected in units of milliliters';
cols.volume_sampled_ml.requirement = 'required'; %
cols.volume_sampled_ml.units       = 'ml'; 
cols.volume_sampled_ml.ecotaxavar  = '';

cols.volume_imaged_ml.description = 'Subset of volume_sampled_ml that was imaged in units of milliliters';
cols.volume_imaged_ml.requirement = 'required'; %
cols.volume_imaged_ml.units       = 'ml'; 
cols.volume_imaged_ml.ecotaxavar  = 'acq_volimage';

cols.pixel_per_um.description = 'Number of pixels per unit length in units of micrometers';    
cols.pixel_per_um.requirement = 'required'; %
cols.pixel_per_um.units       = 'um'; 
cols.pixel_per_um.ecotaxavar  = 'process_pixel'; % Could be process_pixel or acq_pixel 
