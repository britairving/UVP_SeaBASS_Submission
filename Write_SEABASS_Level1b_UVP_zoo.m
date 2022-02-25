function Write_SEABASS_Level1b_UVP_zoo
%% FUNCTION WRITE_SEABASS_LEVEL1B_UVP_ZOO
% Description:
%   Reads in file exported from Ecotaxa Image Module (D.O.I. option) and
%   creates a file for submission to NASA's SEABASS database
%   (https://seabass.gsfc.nasa.gov/wiki/Data_Submission/)
% 
%   The data file of individual particle and plankton level manual
%   classifications is exported from the Ecotaxa Image module in D.O.I.
%   export (tsv) format, and reformatted for submission to SEABASS. The
%   exported tsv file contains variables detailing the acquisition, sample,
%   processing, and ROI. Most variables are in units of pixels and are
%   converted to physical units by the variable process_pixel, the pixel
%   length defined in the calibration document. Depth is corrected for a
%   1.2 m offset due to the difference between the depth sensor and the
%   imaged zone. For the Level 1b data, a ?sample? is defined as each ~1L
%   image acquired during the downcast.  A file containing a list of all
%   samples and their associated depths are included in the associated
%   files in order to maintain provenance, allow the calculation of
%   concentrations, and enable the determination of absence.
%
%   Taxonomic annotation:World Register of Marine Species (WoRMS) is the
%   official taxonomic reference list for the Ocean Biodiversity
%   Information System (OBIS) and was used to provide machine-readable
%   taxonomic identifiers for living organisms. Taxonomic names were
%   matched with a WoRMS AphiaID to the lowest taxonomic rank that could be
%   identified. OBIS can harvest a scientificNameID containing a Life
%   Science Identifier (LSID, http://www.lsid.info), a persistent globally
%   unique identifier for biological objects in a database, from WoRMS. The
%   use of LSIDs from a referenced database permits machine-readable
%   taxonomic identifications. An LSID used in a taxonomic database
%   consists of a uniform resource name (urn) that contains the following
%   information in order: network identifier, a root DNS name of the
%   reference database, a namespace, and the object number or taxon
%   identifier (a living product) that is unique to the biological object
%   defined by that referenced database. An example of a scientificNameID
%   is urn:lsid:marinespecies.org:taxname:345868. In this example,
%   ?urn:lsid? indicates the ID that is specific to life science data and
%   is used for all files, marinespecies.org is the url for the reference
%   database WoRMS, and the namespace ?taxname? informs the user that the
%   following number represents a unique numerical identifier or taxon
%   identifier in WoRMS. In the above example, 345868 represents the taxon
%   identifier (AphiaID) in WoRMS for the Subclass Phaeodaria.
%
%   AphiaID's are used to generate the scientificNameID, and if no AphiaID
%   match is found for a living item, the scientificNameID was generated
%   with the Algae Base id for Eukaryota (i.e.
%   urn:lsid:algaebase.org:taxname:86701).
%
%   Processing scripts publicly available here:  
%   https://github.com/britairving/UVP_submission_formatting 
%
% Steps:
%   1. Read individual particle and plankton level manual
%      classifications exported from the Ecotaxa Image module in D.O.I.
%      export (tsv) format
%   2. Read in metadata from [cruiseid]_UVP_metadata.m
%   3. Read sb fields from SeaBASS_define_taxonomic_Level1b_fields.m
%   4. Query WoRMS database for taxa scientificNames and scientificNameIDs
%   5. Write YAML-formatted namespace file defining non-conforming ROI's
%   6. Calculate sb fields from definitions
%   7. Build SeaBASS formatted metadata header 
%   8. write level1 sb file
%   9. Write Assessed_id_list_[cruiseid].csv file
%
% References:
%  https://github.com/britairving/UVP_submission_formatting
%
%  Picheral M, Colin S, Irisson J-O (2017). EcoTaxa, a tool for the
%  taxonomic classification of images. http://ecotaxa.obs-vlfr.fr.
%
%  OCB_PTWG_TM_FinalDraft_Dec2020.docx
%
%  https://seabass.gsfc.nasa.gov/wiki/plankton_and_particles
%
%  https://github.com/ecotaxa/ecotaxatoolbox/blob/master/UVPread_persample.m
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
%  Brita Irving     <bkirving@alaska.edu>

%% ** USER INPUT REQUIRED ** > Define cruise ID
%cruiseid = 'SR1812';% EXPORTSNP survey cruise R/V Sally Ride 2018
%cruiseid = 'RR1813';% EXPORTSNP process cruise R/V Roger Revelle 2018
cruiseid = 'DY131';% EXPORTSNA survey cruise R/V Discovery 2021

%% ** USER INPUT REQUIRED ** > Define project directory
if ismac 
  projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);
else
  projectdir = fullfile('D:','MATLAB','UVP_project_data',cruiseid);
end

%% Read in basic metadata for project
try
  % Change directory to project folder so can easily call metadata file
  pwd_now = pwd;
  cd(projectdir)
  eval(['hdr = ' cruiseid '_UVP_metadata'])
  % Change directory back to previous working directory
  cd(pwd_now);
catch % catch and explain why script stopped
  fprintf('Cannot load hdr for this project\n')
  fprintf('Need to set up %s_UVP_metadata.m script, see example provided\n',cruiseid)
  keyboard
end
% Test that hdr loaded
if ~exist('hdr','var')
  error('*** Need "hdr" structure - called from [cruiseid]_UVP_metadata.m')
end

%% Level 1b fieldnames
% Review SeaBASS_define_taxonomic_Level1b_fields.m to add new fields or
% adjust fields as necessary. 
fields = SeaBASS_define_taxonomic_Level1b_fields;

%% Configuration
cfg.write_yaml_file = 0; % 1 = writes namespace file (YAML formatted), descriptions in 
cfg.single_sb_file  = 1; % 1 = writes all data to a single sb file. 0 = writes sb file for each eventID (cast in our case)
cfg.ptwg_namespace             = struct();  % Writes associated terms to YAML file
cfg.ptwg_namespace.id          = {'bad_image' 'bead' 'bubble' 'detritus' 'fecal_pellet' 'other'};
cfg.ptwg_namespace.description = 'The Ocean Carbon and Biogeochemistry Phytoplankton Taxonomy Working Group (PTWG) namespace for non-conforming ROIs related to the imaging of plankton and other particles. Version 1.';
cfg.ptwg_namespace.url         = 'https://seabass.gsfc.nasa.gov/ptwg_namespace_v1';

%% Build read filename 
raw_rfile = fullfile(projectdir,hdr.ecotaxaf,[hdr.ecotaxaf '.tsv']); % based on project directory and exported ecotaxa name
assessed_id_file = ['Assessed_id_list_' cruiseid '.csv'];            % Providing a list of all scientificName/scientificNameID pairs assessed by the automated classifier with the data submission enables the determination of both presence and absence of annotations in the Level 1b file. Supplementary lists of which taxonomic categories were assessed by manual and/or automatic classification methods are strongly recommended and are required as part of data submissions if not every ROI in a given datafile was classified. If every ROI was not classified, these lists are essential for the downstream creation of summary products involving the concentrations of phytoplankton taxa. 

%% Unique namespace for non-conforming ROIs
% Define namespace YAML-formatted filename used to define non-conforming
% Use the same namespace for both EXPORTSNP cruises RR1813 and
% SR1812 because submitted at the same time... February 2021
if isfield(hdr,'namespace')
  namespace = hdr.namespace;
else
  namespace = ['namespace_' cruiseid '.yaml'];
end

%% Read raw zoo data
fprintf('\n  Reading raw Zooplankton data from... %s\n',raw_rfile)
% check that file exists
if ~exist(raw_rfile,'file')
  fprintf('File does not exist - make sure path and name are correct: %s\n',raw_rfile)
end
% read raw data from file
raw = readtable(raw_rfile,'FileType','text');
% Convert all doubles to cell array of characters so that
% SeaBASS_define_taxonomic_Level1b_fields.m equations will work
fprintf('Converting numeric columns to strings\n')
for ns = 1:numel(raw.Properties.VariableNames)
  sfield = raw.Properties.VariableNames{ns};
  if isnumeric(raw.(sfield))
    raw.(sfield) = cellstr(num2str(raw.(sfield)));
  end
end

%% Remove bad object_id 
if isfield(hdr,'bad_object_id')
  bad = contains(raw.object_id,hdr.bad_object_id);
  raw(bad,:) = [];
end

%% Pull out unique taxa
% This is necessary to query WoRMS and get the ScientificName and
% ScientficNameID 
[original_fields  ,iu] = unique(raw.object_annotation_hierarchy);
% Initialize fields
raw.child_name  = raw.object_annotation_hierarchy;
raw.parent_name = raw.object_annotation_hierarchy;
for nt = 1:numel(original_fields  )
  str = strsplit(original_fields  {nt},'>');
  idx_nt = strcmp(raw.object_annotation_hierarchy,original_fields  {nt});
  raw.child_name(idx_nt) =  str(end);
  raw.parent_name(idx_nt) = str(end-1);
end

%% Query WoRMS to pull out AphiaID match for each taxonomic name
% save_taxa_filename = [cruiseid 'taxa_' date '.mat'];
save_taxa_filename = fullfile(projectdir,[cruiseid 'taxa.mat']);
if ~exist(save_taxa_filename,'file')
  % Initialize table used as input to WoRMS_AphiaID_taxa_match.m script.
  % Inputs:
  %   taxa | table with fields Name, Name_parent, & Name_original
  %   Example
  %     taxa = table();
  %     taxa.Name          = {};
  %     taxa.Name_parent   = {};
  %     taxa.Name_original = {};
  taxa = table();
  taxa.ID            = {}; % suffix added to abun_, biovol_, avgesd_ fields. This field is not required for WoRMS_AphiaID_taxa_match.m
    
  % Loop through and prase column names new fieldnames with associated ID
  for n_id = 1:numel(original_fields)
    taxa.ID{n_id}   = [num2str(n_id) 'id'];
  end
  taxa.Name          = raw.child_name(iu);                  % taxa name, lowest level.        Required field in WoRMS_AphiaID_taxa_match.m.
  taxa.Name_parent   = raw.parent_name(iu);                 % taxa parent name, if available. Required field in WoRMS_AphiaID_taxa_match.m.
  taxa.Name_original = raw.object_annotation_hierarchy(iu); % name of original variable.      This field is not required for WoRMS_AphiaID_taxa_match.m.

  % Call script to go through and query AphiaID matches
  check_children_manual = false; % true = checks full classifcation of parent if child doesn't have a match and asks user. false = skips this.
  taxa = WoRMS_AphiaID_taxa_match(taxa,check_children_manual,'ecotaxa');
  
  
  % Save to matlab file for simply loading next time
  fprintf('Saving taxa matches to %s\n',save_taxa_filename)
  save(save_taxa_filename,'taxa');
  
else % Load data instead of rerunning everything
  fprintf('Loading..%s\n',save_taxa_filename)
  load(save_taxa_filename);
end

%% Limit data to specific taxa as defined by hdr.limit_to_taxa
% Cell array with taxonomic names selected by the user
% For example:  hdr.limit_to_taxa = {'Salpida' 't008'};   
if isfield(cfg,'limit_to_taxa')
  idx_taxa = ismember(taxa.Name,hdr.limit_to_taxa);
  taxa = taxa(idx_taxa,:);
  idx_raw = ismember(raw.child_name,hdr.limit_to_taxa);
  raw = raw(idx_raw,:);
end

%% Add scientificName and scientificNameID to raw
% initialize fields
raw.scientificName   = raw.child_name; % These will be overwritten, just initializing with correct size & class for now
raw.scientificNameID = raw.child_name; % These will be overwritten, just initializing with correct size & class for now
% Loop through and assign the WoRMS scientificName and scientificNameID
for nt = 1:numel(taxa.Name)
  taxaname = taxa.Name{nt};
  idx_name = find(strcmp(raw.child_name,taxaname));
  % Skip if no matches found
  if isempty(idx_name)
    continue
  end
  % If non-conforming, point to namespace
  if isempty(taxa.scientificName{nt})
    raw.scientificName(idx_name)   = raw.child_name(idx_name);
    raw.scientificNameID(idx_name) = {[namespace ':' raw.child_name{idx_name(1)}]};
  else
    raw.scientificName(idx_name)   = taxa.scientificName(nt);
    raw.scientificNameID(idx_name) = taxa.scientificNameID(nt);
  end
%   % test
%   if strcmp(class(raw.scientificName{idx_name(end)}),'double')
%     fprintf('bad case here... check why this happened\n')
%     keyboard
%   end
end

%% Write YAML-formatted namespace file defining non-conforming ROI's
% These IDs are paired with definitions and are stored in a YAML  formatted
% file in order to serve as a machine-readable configuration file for
% anyone working with the data files.
if cfg.write_yaml_file
  % open file to write
  fid = fopen(fullfile(projectdir,[namespace '.yml']),'w');
  idx_nonconforming = find(contains(raw.scientificNameID,'namespace'));
  [IDs,iu] = unique(raw.scientificNameID(idx_nonconforming));
  iu = idx_nonconforming(iu);
  if any(ismember(raw.scientificName(iu),cfg.ptwg_namespace.id))
    fprintf(fid,'- prefix: ptwg\n');
    fprintf(fid,'  description: %s\n',cfg.ptwg_namespace.description);
    fprintf(fid,'  url: %s\n',cfg.ptwg_namespace.url);
  end
  fprintf(fid,'- prefix: %s\n',namespace);
  fprintf(fid,'  description: \n');
  fprintf(fid,'  url: \n');
  fprintf(fid,'  terms: \n');
  
  for n_id = 1:numel(IDs)
    id = raw.scientificName{iu(n_id)};
    % print to screen formatted text to insert into YAML file
    fprintf(fid,'  - id: %s\n',id);
    if isfield(hdr.nonconforming,id)
      fprintf(fid,'    definition: %s\n',hdr.nonconforming.(id));
    else
      fprintf(fid,'    definition: \n'); % just print
    end
    % See if matches directly to OCB PTWG ids already defined. If there is
    % a match, note the associated term
    if ismember(id,cfg.ptwg_namespace.id)
      fprintf(fid,'    associated_terms: \n');
      fprintf(fid,'    - id: %s\n',['ptwg:' id]);
    end
  end
end

%% Pull in R2R event from Level2 file
% Level2 files are formatted from the detailed ODV file exported from
% Ecotaxa Particle Module. See Write_SEABASS_Level2_UVP_zoo.m and
% associated [cruiseid]_UVP_R2R_merge.m script for specifics.
submitted_files = dir(fullfile(projectdir,'submit','*.sb'));
level2idx = find(contains({submitted_files.name},'TaxonomicLevel2'));
if ~isempty(level2idx)
  level2file = submitted_files(level2idx).name;
  % Read in Level2 file
  % readsb.m from https://seabass.gsfc.nasa.gov/wiki/seabass_tools
  L2 = readsb(fullfile(projectdir,'submit',level2file),'MakeStructure',true);
  if isfield(L2,'r2r_event')
    try %to match raw.sample_id to L2.station_alt_id (cast name equivalent)
      % Initialize r2r_event in raw
      raw.r2r_event = raw.sample_id; % values will be changed
      % pull out unique casts
      ucastname = unique(raw.sample_id);
      for nc = 1:numel(ucastname)
        idx_L2  = find(contains(L2.station_alt_id,ucastname{nc},'IgnoreCase',true));
        if ~isempty(idx_L2)
          r2r_event = L2.r2r_event(idx_L2(1));
          idx_raw = contains(raw.sample_id,ucastname{nc},'IgnoreCase',true);
          raw.r2r_event(idx_raw) = r2r_event;
        end
      end % Loop through unique casts
    catch
      % remove field so does not erroneously write sample_id as r2r_event
      raw.r2r_event = [];
    end
  end % If R2R_EVENT is a field in the Level2 SeaBASS file
end

%% Add all the information from the DAT files in the UVP work folder
% FROM OCB_PTWG_TM_FinalDraft_Dec2020.docx file
% Supplementary lists of which taxonomic categories were assessed by manual
% and/or automatic classification methods are strongly recommended and are
% required as part of data submissions if not every ROI in a given datafile
% was classified. If not every ROI was classified, these lists are
% essential for the downstream creation of summary products involving the
% concentrations of phytoplankton taxa. When every ROI is classified, these
% lists are useful for determining absence. These lists may be specific to
% a given water sample or datafile, e.g., if only diatoms are classified in
% a sample, or they may be comprehensive of every class in a classifier.
if isfield(hdr,'dat_info_file')
  if exist(hdr.dat_info_file,'file')
    % hdr.dat_info_file is a single text file containing all the .dat files
    % merged into a single file from the "work" folders as processed by
    % zooprocess. Generated with UVP_DAT_file_info_merge.m
    % Headers are described in UVP5_User_Manual_20210201 section 2.5.4
    % "DAT files description"
    % Headers in: 
    %   'index' = the image number, starting at 1 in the RAW dat files.
    %   'image' = corresponding image acquisition time (YYYYMMDDHHMMSS_ms)
    %   'p_centibar'  = Pressure in centi-bars
    %   'datfilename' = name of the dat file read in
    %   'work_folder' = name of the folder where dat file was located
    dat = readtable(hdr.dat_info_file,'FileType','text');
    % First, convert pressure centibar to pressure decibar
    dat.p_dbar = dat.p_centibar./10;
    % Match datfilename to raw.sample_profileid to get latitude
    dat.profile_id = lower(erase(dat.datfilename,'.dat'));
    uniqprofile_id = unique(dat.profile_id,'stable');
    % initialize latitude
    dat.latitude  = nan(size(dat.p_dbar));
    % Loop through and get latitude for each profile
    for np = 1:numel(uniqprofile_id)
      pid  = uniqprofile_id(np);
      idx_dat = strcmp(dat.profile_id,pid);
      idx_raw = find(strcmp(raw.sample_profileid,pid),1);
      if isempty(idx_raw)
        continue
      end
      % assign latitude
      if iscell(raw.object_lat)
        dat.latitude(idx_dat) = str2double(raw.object_lat(idx_raw));
      else
        dat.latitude(idx_dat) = raw.object_lat(idx_raw);
      end
    end % UNIQUE PROFILES 
    
    % Calculate depth
    dat.depth = -1.0*gsw_z_from_p(dat.p_dbar,dat.latitude);
    %  Correct depth for the 1.2 m difference in the depth between the two
    %  bases (depth transmitted to ecotaxa images is the depth recorded by
    %  the sensor which is 1.2m above the imaged zone /// ecotaxa particles
    %  does this correction but not the image module)
    dat.depth_corrected =  dat.depth + 1.2;
    
    % %% check difference between depths
    %     % Initialize depth_tsv field - this will be the depth corrected by the
    %     % 1.2m offset.
    %     dat.depth_tsv = nan(size(dat.p_dbar));
    %     % pull out the raw vignette name which should be
    %     % yyyymmddHHMMSS_mmm_xxxx
    %     rawvignettes  = extractBefore(raw.object_rawvig,19);
    %     unique_rawvig = unique(rawvignettes,'stable');
    %     for nimage = 1:numel(unique_rawvig)
    %       rawvig = unique_rawvig(nimage);
    %       idx_raw = find(contains(raw.object_rawvig,rawvig),1);
    %       idx_dat = find(strcmp(dat.image,rawvig));
    %       if ~isempty(idx_dat)
    %         dat.depth_tsv(idx_dat) = raw.object_depth_min(idx_raw)+1.2;
    %       end
    %     end % UNIQUE IMAGES
    %% Write to file
    hdr.merged_dat_file = [cruiseid '_merged_DAT.csv'];
    write_dat_filename = fullfile(projectdir,hdr.ecotaxaf,hdr.merged_dat_file);
    fprintf('Writing data to %s\n', write_dat_filename)
    writetable(dat, write_dat_filename,'FileType','text','Delimiter',',');
  end
end

%% Build data structure, column names, column units, and formatting for file writting
% Variables are either pulled in based on the "rawfield" field, or
% calculated using the eval function, where equations are defined in the
% SeaBASS_define_taxonomic_Level1b_fields.m script.  This was done so
% future users could easily see how variables are calculated.
%
% INITIALIZE
uvp = struct(); % data structure for required/recommended/optional fields
colstrns = {};  % cell array built with names, will be converted to char
colunits = {};  % cell array built with units, will be converted to char
fmt      = {};  % format used to write data to file
% LOOP THROUGH FIELDS DEFINED IN SeaBASS_define_taxonomic_Level1b_fields.m
colnames = fieldnames(fields);
for ncol = 1:numel(colnames)
  fname = colnames{ncol};
  try
    if isfield(fields.(fname),'rawfield')
      uvp.(fname) = raw.(fields.(fname).rawfield);
      colstrns = [colstrns; fname];
      colunits = [colunits; fields.(fname).units];
      if iscell(uvp.(fname))
        fmt = [fmt; '%s'];
      else
        fmt = [fmt; '%.4f'];
      end
    elseif isfield(fields.(fname),'calculate')
      fprintf(' calculating... %s = %s\n',fname,fields.(fname).calculate)
      uvp.(fname) = eval(fields.(fname).calculate);
      colstrns = [colstrns; fname];
      colunits = [colunits; fields.(fname).units];
      if iscell(uvp.(fname))
        fmt = [fmt; '%s'];
      else % set resolution to 4 significant figures
        fmt = [fmt; '%.4f'];
      end
    else
      fprintf('unknown case...%s\n',fname)
    end
  catch
    fprintf(' skipping "%s" for now...\n',fname)
    %keyboard
  end
end

%% Set non-living objects  scientificName_manual to -9999  
uvp.scientificName_manual(contains(raw.object_annotation_hierarchy,'not-living')) = {'-9999'};
%% Remove the hierarchy from data_provider_category_manual
% The field data_provider_category_manual should have the last category.
% For example, not-living>artefact>bubble, should be a bubble.
% This was already done above, secondary method below....
uvp.data_provider_category_manual = raw.child_name;
% % done = 0;
% % while ~done
% %   idx_bad = contains(uvp.data_provider_category_manual,'>');
% %   if ~any(idx_bad) % No more > characters in the field
% %     done = 1;
% %   else 
% %     uvp.data_provider_category_manual(idx_bad) = extractAfter(uvp.data_provider_category_manual(idx_bad),'>');
% %   end
% % end

%% Replace all white space with underscore
% for example: rhizaria like, double sphere, ctenophora sp, rhizaria X sp
uvp_fields = fieldnames(uvp);
for nf = 1:numel(uvp_fields)
  sf = uvp_fields{nf};
  if iscell(uvp.(sf))
    uvp.(sf) = strrep(uvp.(sf),' ','_');
  end
end

%% Now change the folder names in the /images/{folder names} as well
image_dir = fullfile(projectdir,hdr.ecotaxaf,'images');
if ~isfolder(image_dir) && any(contains(raw.child_name,' '))
  fprintf('**ERROR** folders containing images of taxa contains spaces, need to replace these with underscores!\n')
  fprintf('  BUT.... cannot find correct /image folder to automatically change those folder names\n')
  fprintf('  stopped here in keyboard mode .. you need to investigate and build folder path appropriately\n')
  fprintf('  currently looked in fullfile(projectdir,hdr.ecotaxaf,"images") = %s\n',image_dir);
  keyboard
elseif any(contains(raw.child_name,' '))
  image_folders = dir(image_dir);
  image_fnames  = {image_folders.name};
  rename_folders = image_fnames(contains(image_fnames,' '));
  for rn = 1:numel(rename_folders)
    movefile(fullfile(image_dir,rename_folders{rn}),fullfile(image_dir,strrep(rename_folders{rn},' ','_')));
  end
end
  

%% Write single sb file for entire dataset, or split by eventID
if cfg.single_sb_file
  num_sb_files = 1;
else
  eventIDs = unique(uvp.eventID);
  num_sb_files = numel(eventIDs);
  % convert to table, so can index data by eventID
  UVP = struct2table(uvp);
  % remove eventID, and R2R_Event if necessary, from fields
  idx_rm = contains(colstrns,{'eventID' 'R2R_Event'});
  colstrns(idx_rm) = [];
  colunits(idx_rm) = [];
  fmt(idx_rm)      = [];
end

%% Combine column names, units, and format into a single string
colunits = strjoin(colunits,','); 
% Combine column names into a single string
colstr = strjoin(colstrns,','); 
% remove trailing comma
fmt = strjoin(fmt,','); 
% add newline character at the end
fmt = [fmt '\n'];


for nsb = 1:num_sb_files
  if num_sb_files > 1
    eventID_range = strcmp(UVP.eventID,eventIDs(nsb));
    uvp = table2struct(UVP(eventID_range,:),'ToScalar',true);
    eventID_sb = uvp.eventID{1};
    uvp = rmfield(uvp,'eventID');
    if isfield(uvp,'R2R_Event')
      R2R_Event_sb = uvp.R2R_Event{1};
      uvp = rmfield(uvp,'R2R_Event');
    else
      R2R_Event_sb = '';
    end
    add_headers ={['/eventID='  eventID_sb]; ['/R2R_Event='   R2R_Event_sb]};
    
  end
  %% Generate SeaBASS header text
  % Pull out raw filename
  [~,original_file,ext] = fileparts(raw_rfile);
  original_file = [original_file ext];
  
  % pull out max/min latitude, longitude, date, and time
  latmax = num2str(max(uvp.lat));
  latmin = min(uvp.lat);
  lonmax = max(uvp.lon);
  lonmin = min(uvp.lon);
  
  % Find the first and last date/time for SEABASS metadata header
  % since odv2 table has already been sorted by datetime (early on), just
  % pull out first and last
  [~,imin] = min(str2double(uvp.date));
  [~,imax] = max(str2double(uvp.date));
  datemax = uvp.date{imax};
  datemin = uvp.date{imin};
  timemax = uvp.time{imax};
  timemin = uvp.time{imin};
  
  volume_sampled_ml = num2str(str2double(raw.acq_volimage{1})*1000);  % Convert from L to mL
  pixel_per_um      = raw.process_pixel{1};                           % already in micrometer
  
  % Generate seabass filename
  if num_sb_files == 1
    sb_filename = [hdr.raw_wfile '_' datemin '-' datemax '_' hdr.sb_release '.sb'];
  else
    sb_filename = [hdr.raw_wfile '_' datemin '_' erase(timemin,':') '_' hdr.sb_release '.sb'];
  end
 
  hdr_SEABASS={'/begin_header';
    ['/investigators='  hdr.investigators];
    ['/affiliations='   hdr.affiliations];
    ['/contact='        hdr.contact];
    ['/experiment='     hdr.experiment];
    ['/cruise='         hdr.cruise];
    ['/station='        hdr.station];
    '/water_depth=NA';
    ['/data_file_name=' sb_filename];
    ['/documents='      strjoin(hdr.documents.Level1b,',') ',' assessed_id_file];
    ['/data_type='      hdr.data_type];
    ['/data_status='    hdr.data_status.Level1b];
    ['/start_date='     datemin];
    ['/end_date='       datemax];
    ['/start_time='     timemin '[GMT]'];
    ['/end_time='       timemax '[GMT]'];
    ['/north_latitude=' num2str(latmax) '[DEG]'];
    ['/south_latitude=' num2str(latmin) '[DEG]'];
    ['/east_longitude=' num2str(lonmax) '[DEG]'];
    ['/west_longitude=' num2str(lonmin) '[DEG]'];
    ['/missing='        hdr.missing];
    ['/delimiter='      hdr.delimiter];
    ['/instrument_manufacturer=' hdr.inst_mfr];
    ['/instrument_model='        hdr.inst_model];
    ['/calibration_files='       hdr.calfiles];
    ['/calibration_date='        hdr.caldates];...
    ['/associated_files='  strcat(hdr.raw_wfile,'_images.tar.gz') ',' original_file];...
    '/associated_file_types=planktonic,raw';...
    ['/volume_sampled_ml=' volume_sampled_ml];...
    ['/volume_imaged_ml='  volume_sampled_ml];...
    ['/pixel_per_um='      pixel_per_um];...
    ['/associatedMedia_source=' hdr.ecotaxa_url]; % e.g. /associatedMedia_source=https://ecotaxa.obs-vlfr.fr/prj/1591
    '/length_representation_instrument_varname=object_major_MULTIPLYBY_process_pixel';... % length_representation_instrument_varname (um): the instrument?s variable name equivalent to ?length_representation? (e.g., maxFeretDiameter).
    '/width_representation_instrument_varname=object_minor_MULTIPLYBY_process_pixel';...  % width_representation_instrument_varname (um): the instrument?s variable name equivalent to ?width_representation? (e.g., minFeretDiameter).
    '!'};
  % Add eventID and R2R_Event to header 
  if num_sb_files > 1
    idx_add_after = find(contains(hdr_SEABASS,'west_longitude'));
    hdr_SEABASS = [hdr_SEABASS(1:idx_add_after); add_headers; hdr_SEABASS(idx_add_after+1:end)];
  end
  % Insert comments, then finish with /fields and /units
  hdr_SEABASS = [hdr_SEABASS; hdr.comments.Level1b;...
    '!';...
    ['/fields=' colstr];...
    ['/units='  colunits];...
    '/end_header'];
  
  % Add associated file with merged DAT information.
  if isfield(hdr,'merged_dat_file')
    hdr_SEABASS(contains(hdr_SEABASS,'associated_files='))      = strcat(hdr_SEABASS(contains(hdr_SEABASS,'associated_files=')),',',hdr.merged_dat_file);
    hdr_SEABASS(contains(hdr_SEABASS,'associated_file_types=')) = strcat(hdr_SEABASS(contains(hdr_SEABASS,'associated_file_types=')),',metadata');
  end
  
  % check if there is whitespace in any metadata headers
  % whitespace in comments is OKAY
  if any(contains(hdr_SEABASS,' ') & ~contains(hdr_SEABASS,'!'))
    fprintf('White space in metadata header, must remove to pass fcheck\n')
    keyboard
  end

  %% Format odv2 table to temporary variable to enable simply write to file
  raw_write = table2cell(struct2table(uvp)); % convert to cell array to handle different variable types
  raw_write = raw_write';            % transpose because fprintf function prints data columnwise
  % convert NaNs to missing identifier
  % raw_write(cellfun(@(x) isnumeric(x) && isnan(x), raw_write)) = {''}; % convert NaNs to '' for raw format
  raw_write(cellfun(@(x) isnumeric(x) && isnan(x), raw_write)) = {str2double(hdr.missing)}; % missing=-9999 SeaBASS required header field
  
  %% Write raw table to file
  if exist(sb_filename,'file')
    fprintf('\n  Deleting data file: %s\n',sb_filename);
    delete(sb_filename)
  end
  
  fprintf('\n  Writing data in raw format to: %s\n',sb_filename);
  fileID = fopen(sb_filename,'w');  % open file
  if fileID < 0
    fprintf(' *** Error opening file %s\n',sb_filename)
    keyboard
  end
  fprintf(fileID,'%s\n',hdr_SEABASS{:});% write header
  fprintf(fileID,fmt,raw_write{:});     % write data
  fclose(fileID);                       % close file
  
 %% Write a smaller version for testing with fcheck
 if num_sb_files == 1
   raw_write_small = raw_write(:,1:200:end);
   fprintf('\n  Writing data in raw format to: %s\n',strrep(sb_filename,'.sb','_fcheck.sb'));
   fileID = fopen(strrep(sb_filename,'.sb','_fcheck.sb'),'w');  % open file
   if fileID < 0
     fprintf(' *** Error opening file %s\n',strrep(sb_filename,'.sb','_fcheck.sb'))
     keyboard
   end
   fprintf(fileID,'%s\n',hdr_SEABASS{:});% write header
   fprintf(fileID,fmt,raw_write_small{:});     % write data
   fclose(fileID);                       % close file
   % To troubleshoot why not printing correctly, comment above and just print to screen
   % fprintf('%s\n',hdr_SEABASS{:}) % write header
   % fprintf(fmt,raw_write{:})      % write data
 end
end
%% List of assessed IDs for automated and/or manual classification 
% Providing a list of all scientificName/scientificNameID pairs assessed by
% the automated classifier with the data submission enables the
% determination of both presence and absence of annotations in the Level 1b
% file. Supplementary lists of which taxonomic categories were assessed by
% manual and/or automatic classification methods are strongly recommended
% and are required as part of data submissions if not every ROI in a given
% datafile was classified. If every ROI was not classified, these lists are
% essential for the downstream creation of summary products involving the
% concentrations of phytoplankton taxa.

assessed = struct();
[~,idx_unique] = unique(uvp.data_provider_category_manual,'stable');
assessed.data_provider_category_manual = uvp.data_provider_category_manual(idx_unique);
if isfield(uvp,'data_provider_category_automated')
  assessed.data_provider_category_automated = uvp.data_provider_category_automated(idx_unique);
end
assessed.scientificName_manual   = uvp.scientificName_manual(idx_unique);
assessed.scientificNameID_manual = uvp.scientificNameID_manual(idx_unique);
% add hierarchy to assess file because it was not allowed to stay in the
% "data_provider_category_manual" field.
assessed.data_provider_category_manual_annotation_hierarchy = raw.object_annotation_hierarchy(idx_unique);
% assessed.data_provider_category_manual_classification = 
assessed = struct2table(assessed);
writetable(assessed,fullfile(projectdir,assessed_id_file),'Delimiter',',');

%% END OF MAIN FUNCTION - SUBFUNCTIONS BELOW
fprintf('workflow is finished\n')

end %% MAIN FUNCTION








