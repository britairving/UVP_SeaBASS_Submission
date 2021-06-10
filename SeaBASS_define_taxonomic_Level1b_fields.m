function fields = SeaBASS_define_taxonomic_Level1b_fields
%% FUNCTION SEABASS_DEFINE_TAXONOMIC_LEVEL1B_FIELDS
% Description:
%   Defines fieldnames for Level1b SeaBASS files. Data are not actually
%   calculated in this script, just defined.
%
%   Each variable has either "rawfield" or "calculate" as a field. This was
%   done so future users could easily see how variables are calculated and
%   edit and/or add fields as necessary.
%     "rawfield"  | field/column name from TSV file exported from Ecotaxa
%     "calculate" | formula for calculation
%
%   For example: variable "lat"
%      fields.lat.rawfield  = 'object_lat'; 
%             would lead to | lat = raw.object_lat; 
%      fields.lat.calculate = "str2double(raw.object_lat)";  
%             would lead to | lat = eval(fields.lat.calculate);
%
% References:
%
%  OCB_PTWG_TM_FinalDraft_Dec2020.docx
%
%  SeaBASS Plankton and other Particles (IFCB, UVP)
%   <https://seabass.gsfc.nasa.gov/wiki/plankton_and_particles>
%
%  SeaBASS Standardized Fields and Units
%   <https://seabass.gsfc.nasa.gov/wiki/stdfields>
%
%  Ecotaxa field definitions
%   <https://sites.google.com/view/piqv/ecotaxa?authuser=0>
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
fields = struct();

fields.station.description = 'Sample station';
fields.station.requirement = 'optional'; %
fields.station.units       = 'none'; 
fields.station.rawfield    = 'sample_stationid'; % e.g. "ss1" "ls_18" "wire_walker_calibration_-_3"

fields.station_alt_id.description = 'Alternate sample station identifier (use ONLY if station has already been used as a field name)';
fields.station_alt_id.requirement = 'optional';
fields.station_alt_id.units       = 'none';
fields.station_alt_id.rawfield    = 'sample_id'; % e.g. ctd001


fields.eventID.description = 'A unique identifier associated with the sample as an event';
fields.eventID.requirement = 'optional'; %
fields.eventID.units       = 'none'; 
fields.eventID.rawfield    = 'object_rawvig'; % This is the same as station_alt_id in the Level2 ZOO file


% Fields not listed in SeaBASS wiki on plankton_and_particles (see
% comments) but important to note
fields.depth.description = 'Depth of measurement';
fields.depth.requirement = 'optional'; %
fields.depth.units       = 'm'; 
fields.depth.calculate   = 'raw.object_depth_min+1.2'; % Note that there is a 1.2 m difference in the depth between the two bases (depth transmitted to ecotaxa images is the depth recorded by the sensor which is 1.2m above the imaged zone /// ecotaxa particles does this correction but not the image module)

fields.lat.description = 'Sample latitude (decimal fractions; -90 to 90)';
fields.lat.requirement = 'optional'; %
fields.lat.units       = 'degrees'; 
fields.lat.calculate   = "str2double(raw.object_lat)"; %

fields.lon.description = 'Sample longitude (decimal fractions; -180 to 180)';
fields.lon.requirement = 'optional'; %
fields.lon.units       = 'degrees'; 
fields.lon.calculate   = "str2double(raw.object_lon)"; % E.g. lon = eval(fields.lon.calculate);

fields.date.description = 'Sample date';
fields.date.requirement = 'optional'; 
fields.date.units       = 'yyyymmdd'; 
fields.date.rawfield    = 'object_date'; % E.g. date = raw.object_date;

fields.time.description = 'Sample time';
fields.time.requirement = 'optional'; 
fields.time.units       = 'hh:mm:ss'; 
fields.time.calculate   = "insertAfter(insertAfter(extractBetween(raw.object_rawvig,9,14),2,':'),5,':')"; % pull out time from object_rawvig field then convert from hhmmss to hh:mm:ss 
% fields.time.calculate   = "insertAfter(insertAfter(raw.object_time,2,':'),5,':')"; % Convert from hhmmss to hh:mm:ss

fields.R2R_Event.description = 'Rolling Deck to Repository (R2R) Unique sample event number';
fields.R2R_Event.requirement = 'optional'; 
fields.R2R_Event.units       = 'none'; 
fields.R2R_Event.rawfield    = 'r2r_event';


fields.associatedMedia.description = 'A unique persistent identifier of the media associated with the occurrence. The field provides the unique imagery file name corresponding to the source of the ROI. Alternately, use this field to provide a URL pointing to a permanent landing page for the ROI image. In the latter case, instructions should be provided as comments in the header on how to construct the local file name based on the URL. If the local imagery file name cannot be constructed from the URL, then list both filenames. Use the pipe character ?|? to separate the names, and do not use spaces.';
fields.associatedMedia.requirement = 'required';
fields.associatedMedia.units       = 'none';
fields.associatedMedia.rawfield    = 'img_file_name'; % e.g. images/Annelida/106749790_0.jpg - path to the image

fields.associatedMedia_source.description = 'a unique persistent URL pointing to the landing page for a water sample from which multiple ROIs are derived';
fields.associatedMedia_source.requirement = 'optional';
fields.associatedMedia_source.units       = 'none';
fields.associatedMedia_source.calculate   = 'extractBefore(raw.object_rawvig,19)'; % e.g. '20180811192720_032' from '20180811192720_032_0001'

% From Dr. Aimee Neeley
% I think Brita has the right idea generally. I would suggest she put
% Aulatractus in the data_provider_category_manual and then Aulosphaeridae
% the ScientificName and ScientificNameID. The ScientificName and
% ScientificNameID should match which is why I would not suggest putting
% Aulatractus as a ScientificName when there isn?t an appropriate AphiaID
% and, therefore, an appropriate LSID.
% I would label them this way:
% data_provider_category_manual = 'Aulatractus'
% scientificName_automated      = 'Aulosphaeridae'
% scientificNameID_automated    = 'urn:lsid:marinespecies.org:taxname:367360'
% scientificName_manual         = 'Aulosphaeridae'
% scientificNameID_manual       = 'urn:lsid:marinespecies.org:taxname:367360'
fields.data_provider_category_manual.description = 'A category used by the data provider to name the organism or particle for a manual identification, not necessarily a scientific name.';
fields.data_provider_category_manual.requirement = 'optional'; % (recommended but optional) 
fields.data_provider_category_manual.units       = 'none';
fields.data_provider_category_manual.rawfield    = 'object_annotation_category';  % e.g. 'Annelida' from object_annotation_hierarchy = 'living>Eukaryota>Opisthokonta>Holozoa>Metazoa>Annelida'

fields.data_provider_category_automated.description = 'A category used by the data provider to name the organism or particle for an automated classification, not necessarily a scientific name (e.g., pennate or detritus).';
fields.data_provider_category_automated.requirement = 'optional'; % (recommended but optional) 
fields.data_provider_category_automated.units       = 'none';
fields.data_provider_category_automated.rawfield    = 'object_annotation_hierarchy'; % e.g. 'living>Eukaryota>Opisthokonta>Holozoa>Metazoa>Annelida'

% The scientificName/scientificNameID pairs can be determined manually by
% searching WoRMS or automatically using web services with a script or with
% the WoRMS Taxon Match Graphical User Interface (GUI) .
fields.scientificName_automated.description = 'A scientific name from a recognized taxonomic reference database (e.g., WoRMS, AlgaeBase) at the lowest level that matches the data providers category for an automated classification paired to a scientificNameID. Generally, the ROI corresponds to an occurrence assigned to a single taxonomic name.'; 
fields.scientificName_automated.requirement = 'optional'; % (recommended but optional) 
fields.scientificName_automated.units       = 'none';
fields.scientificName_automated.rawfield    = 'scientificName';     % Not sure about this....

fields.scientificNameID_automated.description = 'A life science identifier (LSID) from a recognized taxonomic reference database (e.g., WoRMS, AlgaeBase) at the lowest level that matches the data provider?s category for an automated classification. e.g., urn:lsid:marinespecies.org:taxname:233015 where "urn:lsid" indicates the ID that is specific to life science data and is used for all files, marinespecies.org is the url for the reference database WoRMS, and the namespace ?taxname? informs the user that the following number represents a unique numerical identifier or taxon identifier in WoRMS. 233015 represents the taxon identifier (AphiaID) in WoRMS for the dinoflagellate species Karenia brevis.';
fields.scientificNameID_automated.requirement = 'required'; 
fields.scientificNameID_automated.units       = 'none';
fields.scientificNameID_automated.rawfield    = 'scientificNameID'; % Not sure about this....

% fields.scientificName_manual.description = 'A scientific name from a recognized taxonomic reference database (e.g., World Register of Marine Species, AlgaeBase) at the lowest level that matches the data providers category, for a manual identification matched to scientificNameID. Generally, the ROI corresponds to an occurrence assigned to a single taxonomic name.';
% fields.scientificName_manual.requirement = 'optional'; % (recommended but optional) 
% fields.scientificName_manual.units       = 'none';
% fields.scientificName_manual.rawfield    = 'scientificName';
% 
% fields.scientificNameID_manual.description = 'An LSID from a recognized taxonomic reference database (e.g., World Register of Marine Species, AlgaeBase) at the lowest level that matches the data providers category for a manual identification.';
% fields.scientificNameID_manual.requirement = 'required'; %
% fields.scientificNameID_manual.units       = 'none';
% fields.scientificNameID_manual.rawfield    = 'scientificNameID';

fields.area_cross_section.description = 'Cross-sectional area of the target detected within the ROI determined by means specified in the image processing method or protocol document.';
fields.area_cross_section.requirement = 'required'; %
fields.area_cross_section.units       = 'um^2';  
fields.area_cross_section.calculate   = "raw.object_area .* str2double(raw.process_pixel)";

fields.length_representation.description = 'Representation of length of the target detected within the ROI or largest mesh size for which the target could be retained, determined by means specified in the image processing method or protocol document.';
fields.length_representation.requirement = 'required'; %
fields.length_representation.units       = 'um'; 
fields.length_representation.calculate   = "raw.object_major .* str2double(raw.process_pixel)"; % object_major(primary axis of the best fitting ellipse for the object in pixel)*process_pixel(pixel size in ?m)

fields.width_representation.description = 'Representation of width of the target detected within the ROI or smallest mesh size through which the target could pass, determined by means specified in the image processing method or protocol document.';
fields.width_representation.requirement = 'required'; %
fields.width_representation.units       = 'um'; 
fields.width_representation.calculate   = "raw.object_minor .* str2double(raw.process_pixel)";  % object_minor(secondary axis of the best fitting ellipse for the object in pixel)*process_pixel(pixel size in ?m)

fields.equivalent_spherical_diameter.description = 'Equivalent spherical diameter of the target detected within the ROI determined by means specified in the image processing method or protocol document.';
fields.equivalent_spherical_diameter.requirement = 'optional'; %
fields.equivalent_spherical_diameter.units       = 'um'; 
fields.equivalent_spherical_diameter.calculate   = "raw.object_esd .* str2double(raw.process_pixel)"; % object_esd(!!!assuming this is equivalent spherical diameter in pixels!!!)*process_pixel(pixel size in ?m)

% fields.area_based_diameter.description = 'Area-based diameter of the target detected within the ROI determined by means specified in the image processing method or protocol document.';      
% fields.area_based_diameter.requirement = 'optional'; %
% fields.area_based_diameter.units       = 'um'; 
% fields.area_based_diameter.calculate   = "2 .* sqrt(uvp.area_cross_section ./ pi)"; % double the radius, where radius = sqrt(area/pi)

fields.biovolume.description = 'Biovolume for the target detected within the ROI determined by means specified in the biovolume calculation method or protocol document.';
fields.biovolume.requirement = 'required'; %
fields.biovolume.units       = 'um^3';  
fields.biovolume.calculate   = '(pi/6).*uvp.length_representation.*uvp.width_representation.^2'; % Let's calculate biovolume for each ROI by using major and minor axis dimensions from Ecotaxa (object_major, object_minor) and assuming they are prolate spheroids.  In that case, biovolume=(pi/6)*object_major*object_minor^2.

% fields.biovolume.calculate   = "4/3 .* pi .* uvp.area_based_diameter.^3";           % Need to figure this out!!!
% fields.biovolume.calculate   = "raw.object_circ_.^3 ./ (6*pi^2)";
% fields.biovolume.calculate   = "4/3 .* pi .* uvp.equivalent_spherical_diameter.^3"; % Need to figure this out

% Below are Required in the file header or as a field
% For the UVP put volume_imaged_ml and volume_sampled_ml in header because
% they are the same
% fields.volume_imaged_ml.description = 'Subset of volume_sampled_ml that was imaged in units of milliliters';
% fields.volume_imaged_ml.requirement = 'required'; %
% fields.volume_imaged_ml.units       = 'ml'; 
% fields.volume_imaged_ml.calculate   = "uvp.volume_sampled_ml"; 

% % These should always be the same, so just put in metadata headers
% fields.acq_volimage.description = 'Original volume of sample collected in units of milliliters';
% fields.volume_sampled_ml.requirement = 'required'; %
% fields.volume_sampled_ml.units       = 'ml'; 
% fields.volume_sampled_ml.calculate   = "str2double(raw.acq_volimage) ./ 1000.0"; % Convert from L to mL
% 
% fields.pixel_per_um.description = 'Number of pixels per unit length in units of micrometers';    
% fields.pixel_per_um.requirement = 'required'; %
% fields.pixel_per_um.units       = 'um'; 
% fields.pixel_per_um.calculate    = "str2double(raw.process_pixel)"; % Could be process_pixel or acq_pixel 

