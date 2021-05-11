function Write_SEABASS_Level1b_UVP_zoo
%% FUNCTION WRITE_SEABASS_LEVEL1B_UVP_ZOO
% Description:
%   Reads in detailed raw formatted file exported from Ecotaxa Particle
%   Module and creates a file for submission to NASA's SEABASS database
%   https://seabass.gsfc.nasa.gov/wiki/Data_Submission/
%
% Steps:
%   1. Read in detailed raw PAR file
%   2. Format fields and column names following SeaBASS format
%   3. Define metadata header that will be printed to .sb file
%   4. Write sb file with DNSD and DVSD
%   5. Write suplementary sb file with non-differential NSD and VSD
%
% References:
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
%  Brita Irving <bkirving@alaska.edu>
%  Andrew McDonnell <amcdonnell@alaska.edu>
%% Configuration
cfg.write_yaml_file = 0;                    % 1 = writes namespace file (YAML formatted), descriptions in 
cfg.limit_to_taxa   = {'Salpida' 't008'};   % Cell array with limited taxa names
cfg.ptwg_namespace             = struct();  % Writes associated terms to YAML file
cfg.ptwg_namespace.id          = {'bad_image' 'bead' 'bubble' 'detritus' 'fecal_pellet' 'other'};
cfg.ptwg_namespace.description = 'Ocean Carbon and Biogeochemistry Phytoplankton Taxonomy Working Group';
cfg.ptwg_namespace.url         = 'https://seabass.gsfc.nasa.gov/wiki/ptwg_namespace';
%% Level 1b fieldnames
fields = SeaBASS_define_taxonomic_Level1b_fields;

fprintf('LOOK THROUGH https://github.com/ecotaxa/ecotaxatoolbox/blob/master/UVPread_persample.m\n')
fprintf('CALCULATE NECESSARY VARIABLES...\n');

%% Define read and write filenames
%% EXPORTSNP survey cruise R/V Sally Ride
cruiseid = 'SR1812';
ecotaxaf = 'ecotaxa_export_1591_20210329_1704'; % 'ecotaxa_export_1591_20210119_0944';%'task_42058_export_1591_20210329_1704'; % Folder name of ecotaxa exported RAW data
projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);

raw_wfile = 'EXPORTS-EXPORTSNP_UVP5-TaxonomicLevel1b_Salps_survey_20180814-20180909_R0.sb';

%% EXPORTSNP process cruise R/V Roger Revelle
% cruiseid = 'RR1813';
% projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);
% raw_wfile = 'EXPORTS-EXPORTSNP_UVP5-TaxonomicLevel1b_process_20180814-20180909_R0.sb';
% r2r_elog  = fullfile(projectdir,'R2R_ELOG_RR1813_FINAL_EVENTLOG_20180912_170812.xlsx');

%% Build read filename based on project directory and exported ecotaxa name
raw_rfile = fullfile(projectdir,ecotaxaf,[ecotaxaf '.tsv']);

%% Unique namespace for non-conforming ROIs
% Define namespace YAML-formatted filename used to define non-conforming
% Use the same namespace for both EXPORTSNP cruises RR1813 and
% SR1812 because submitted at the same time... February 2021
switch cruiseid
  case {'SR1812' 'RR1813'}
    namespace = 'namespace_EXPORTSNP_UVP';
  otherwise
    namespace = ['namespace_' cruiseid '.yaml'];
end

%% Specify write filename
% Release number R0 suggested if data_type=preliminary
 % specify write file, this is not computer specific because uses same path
% Write files to submit folder
submit_dir = fullfile(projectdir,'submit');
if ~isdir(submit_dir)
  try
    mkdir(submit_dir);
  catch
    error('cannot create "submit" directory, maybe a permissions error?')
  end
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

%% Read raw zoo data
fprintf('\n  Reading raw Zooplankton data from... %s\n',raw_rfile)

% check that file exists
if ~exist(raw_rfile,'file')
  fprintf('File does not exist - make sure path and name are correct: %s\n',raw_rfile)
end
% read raw data from file
raw = readtable(raw_rfile,'FileType','text');


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
  fprintf('Saving taxa matches to %s\n',save_taxa_filename)
  save(save_taxa_filename,'taxa');
else % Load data instead of rerunning everything
  fprintf('Loading..%s\n',save_taxa_filename)
  load(save_taxa_filename);
end

%% Limit data to specific taxa as defined by cfg.limit_to_taxa
% Cell array with taxonomic names selected by the user
% For example:  cfg.limit_to_taxa = {'Salpida' 't008'};   
if isfield(cfg,'limit_to_taxa')
  idx_taxa = ismember(taxa.Name,cfg.limit_to_taxa);
  taxa = taxa(idx_taxa,:);
  idx_raw = ismember(raw.child_name,cfg.limit_to_taxa);
  raw = raw(idx_raw,:);
end

%% Add scientificName and scientificNameID to raw
% initialize fields
raw.scientificName   = raw.child_name; % These will be overwritten
raw.scientificNameID = raw.child_name; % These will be overwritten
for nt = 1:numel(taxa.Name)
  taxaname = taxa.Name{nt};
  idx_name = find(strcmp(raw.child_name,taxaname));
  raw.scientificName(idx_name)   = taxa.scientificName(nt);
  raw.scientificNameID(idx_name) = taxa.scientificNameID(nt);
  % If non-conforming, point to namespace
  if isempty(taxa.scientificName{nt})
    raw.scientificName(idx_name)   = raw.child_name(idx_name);
    raw.scientificNameID(idx_name) = {[namespace ':' raw.child_name{idx_name(1)}]};
  end
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
    fprintf(fid,'- prefix: ptwg\n')
    fprintf(fid,'  description: %s\n',cfg.ptwg_namespace.description);
    fprintf(fid,'  url: %s\n',cfg.ptwg_namespace.url);
  end
  fprintf(fid,'- prefix: %s\n',namespace)
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
submitted_files = dir(fullfile(submit_dir,'*.sb'));
level2idx = find(contains({submitted_files.name},'TaxonomicLevel2'));
if ~isempty(level2idx)
  level2file = submitted_files(level2idx).name;
  % Read in Level2 file
  % readsb.m from https://seabass.gsfc.nasa.gov/wiki/seabass_tools
  L2 = readsb(fullfile(submit_dir,level2file),'MakeStructure',true);
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
fmt      = '';  % format used to write data to file
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
        fmt = [fmt '%s,'];
      else
        fmt = [fmt '%f,'];
      end
    elseif isfield(fields.(fname),'calculate')
      fprintf(' calculating... %s = %s\n',fname,fields.(fname).calculate)
      uvp.(fname) = eval(fields.(fname).calculate);
      colstrns = [colstrns; fname];
      colunits = [colunits; fields.(fname).units];
      if iscell(uvp.(fname))
        fmt = [fmt '%s,'];
      else
        fmt = [fmt '%f,'];
      end
    else
      fprintf('unknown case...%s\n',fname)
    end
  catch
    fprintf(' skipping "%s" for now...\n',fname)
    %keyboard
  end
end

% Combine all units into a single string
colunits = strjoin(colunits,','); 
% Combine column names into a single string
colstr = strjoin(colstrns,','); 
% remove trailing comma
fmt(end) = [];
% add newline character at the end
fmt = [fmt '\n'];


%% Build filename for assessed IDs
assessed_id_file = ['Assessed_id_list_' cruiseid '.csv'];

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

volume_sampled_ml = num2str(str2double(raw.acq_volimage{1})/1000);  % Convert from L to mL
pixel_per_um      = raw.process_pixel{1};                           % already in micrometer

hdr_SEABASS={'/begin_header';
  ['/investigators='  hdr.investigators];
  ['/affiliations='   hdr.affiliations];
  ['/contact='        hdr.contact];
  ['/experiment='     hdr.experiment];
  ['/cruise='         hdr.cruise];
  ['/station='        hdr.station];
  ['/data_file_name=' raw_wfile];
  ['/associated_files=' original_file]; % The original name of the data file exported from Ecotaxa
  '/associated_file_types=raw';
  ['/documents='      strjoin(hdr.documents.Level1b,',')];
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
  ['/instrument_manufacturer='  hdr.inst_mfr];
  ['/instrument_model='         hdr.inst_model];
  ['/calibration_files='        hdr.calfiles];
  ['/calibration_date='         hdr.caldates];...
  ['/associated_files=images,'  assessed_id_file ',' original_file];...
  '/associated_file_types=planktonic,Assessed_IDs_list,raw';...
  ['/volume_sampled_ml=' volume_sampled_ml];...
  ['/pixel_per_um='      pixel_per_um];...
  'length_representation_instrument_varname=object_major*process_pixel';... % length_representation_instrument_varname (um): the instrument?s variable name equivalent to ?length_representation? (e.g., maxFeretDiameter).  
  'width_representation_instrument_varname=object_minor*process_pixel';...  % width_representation_instrument_varname (um): the instrument?s variable name equivalent to ?width_representation? (e.g., minFeretDiameter).   
  '!'};

% Insert comments, then finish with /fields and /units
hdr_SEABASS = [hdr_SEABASS; hdr.comments.Level1b;...
  '!';
  ['/fields=' colstr];
  ['/units='  colunits];
  '/end_header'];

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
if exist(raw_wfile,'file')
  fprintf('\n  Deleting data file: %s\n',raw_wfile);
  delete(raw_wfile)
end

fprintf('\n  Writing data in raw format to: %s\n',fullfile(submit_dir,raw_wfile));
fileID = fopen(fullfile(submit_dir,raw_wfile),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile(submit_dir,raw_wfile))
  keyboard
end
fprintf(fileID,'%s\n',hdr_SEABASS{:});% write header
fprintf(fileID,fmt,raw_write{:});     % write data
fclose(fileID);                       % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
fprintf('%s\n',hdr_SEABASS{:}) % write header
fprintf(fmt,raw_write{:})      % write data

%% List of assessed IDs for automated and/or manual classification 
assessed = struct();
[~,idx_unique] = unique(uvp.data_provider_category_manual,'stable');
assessed.data_provider_category_manual    = uvp.data_provider_category_manual(idx_unique);
assessed.data_provider_category_automated = uvp.data_provider_category_automated(idx_unique);
assessed.scientificName_automated   = uvp.scientificName_automated(idx_unique);
assessed.scientificNameID_automated = uvp.scientificNameID_automated(idx_unique);

assessed = struct2table(assessed);
writetable(assessed,fullfile(projectdir,assessed_id_file),'Delimiter',',');

%% END OF MAIN FUNCTION - SUBFUNCTIONS BELOW
fprintf('workflow is finished\n')

end %% MAIN FUNCTION








