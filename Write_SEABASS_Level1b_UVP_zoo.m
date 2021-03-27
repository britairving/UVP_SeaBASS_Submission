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
%% Level 1b fieldnames
cols = SeaBASS_define_taxonomic_Level1b_fields;

fprintf('LOOK THROUGH https://github.com/ecotaxa/ecotaxatoolbox/blob/master/UVPread_persample.m\n')
fprintf('CALCULATE NECESSARY VARIABLES...\n');
keyboard
%% Define read and write filenames
%% EXPORTSNP survey cruise R/V Sally Ride
cruiseid = 'SR1812';
ecotaxaf = 'ecotaxa_export_1591_20210119_0944'; % Folder name of ecotaxa exported RAW data
projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);

raw_wfile = 'EXPORTS-EXPORTSNP_UVP5-TaxonomicLevel1b_survey_20180814-20180909_R0.sb';

%% EXPORTSNP process cruise R/V Roger Revelle
% cruiseid = 'RR1813';
% projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);
% raw_wfile = 'EXPORTS-EXPORTSNP_UVP5-TaxonomicLevel1b_process_20180814-20180909_R0.sb';
% r2r_elog  = fullfile(projectdir,'R2R_ELOG_RR1813_FINAL_EVENTLOG_20180912_170812.xlsx');

%% Build read filename based on project directory and exported ecotaxa name
raw_rfile = fullfile(projectdir,ecotaxaf,[ecotaxaf '.tsv']);
namespace = ['namespace_' cruiseid '.yaml'];
%%
% Define namespace YAML-formatted filename used to define non-conforming
% Use the same namespace for both EXPORTSNP cruises RR1813 and
% SR1812 because submitted at the same time... February 2021
switch cruiseid
%   case {'SR1812' 'RR1813'}
%     namespace = 'EXPORTSNP_Ecotaxa';
  otherwise
    namespace = ['EXPORTS_' cruiseid];
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

fields = fieldnames(raw);
for f = 1:numel(fields)
  fprintf('%s\n',fields{f});
end
keyboard
%% Remove unwanted columns
% define variables that we want to remove
remove_fields = {'DataOwner' 'Cruise' 'Rawfilename' 'Instrument' 'SN' 'CTDrosette'};
remthis       = find(contains(cols2,remove_fields));
odv2(:,remthis) = [];
cols2(remthis)  = [];

%% Create date and time variables and sort by datetime
fnames = fieldnames(odv2);
% Convert datetime variable to string
idx_date_field = find(contains(cols2,'yyyy'));
date_field = fnames{idx_date_field};
if isdatetime(odv2.(date_field))
  % sort odv2 table by datetime field
  odv2      = sortrows(odv2,idx_date_field);
  odv2.date = datestr(odv2.(date_field),'yyyymmdd');
  odv2.time = datestr(odv2.(date_field),'HH:MM:SS');
  % now convert original date_field to string
  odv2.(date_field) = datestr(odv2.(date_field),'yyyy-mm-dd HH:MM:SS');
  cols2{1,end+1}='date';
  cols2{1,end+1}='time';
else % add catch so can fix before errors are generated to avoid confusion
  fprintf('Unexpected format! Stopping here so can figure out why datetime variable is not in datetime format\n')
  keyboard
end

%% switch order of columns and remove datetime
ncols = numel(odv2.Properties.VariableNames);
% find index location of variables we want to reorder
i_datetime = find(strcmp(odv2.Properties.VariableNames,'datetime'));
i_date     = find(strcmp(odv2.Properties.VariableNames,'date'));
i_time     = find(strcmp(odv2.Properties.VariableNames,'time'));
new_order  = [1:i_datetime i_date i_time i_datetime+1:ncols];
% Reorder cell array with full column names
cols2 = cols2(new_order); % replaced cols2=[cols2(1,1:dt-1) , cols2(1,d) , cols2(1,d+1) , cols2(1,dt+1:d-1) , cols2(1,d+2:end)];
% added catch for older versions of MATLAB
try % movevars introduced in R2018a
  odv2=movevars(odv2,'date','After','datetime');
  odv2=movevars(odv2,'time','After','date');
catch % if older version of MATLAB
  % this will work for all versions of matlab but kept the above method for
  % clarity
  odv2 = odv2(:,new_order);
end
% sort by datetime then remove the datetime column
odv2.datetime     = [];
cols2(i_datetime) = [];
% delete redundant variables that result from reordering
odv2.date_1 = [];
odv2.time_1 = [];
cols2(end-1:end) = []; % date and time at the end


%% Query WoRMS to pull out AphiaID match for each taxonomic name
% save_taxa_filename = [cruiseid 'taxa_' date '.mat'];
save_taxa_filename = [cruiseid 'taxa.mat'];
if ~exist(save_taxa_filename,'file')
  % First step | parse taxonomic names from Ecotaxa fields
  % Only need to look at abundance names, because column names are the same
  % for biovolume and avgesd.
  original_fields = cols2(contains(cols2,'[#/m^3]'));
  original_fields = erase(original_fields,'_[#/m^3]');
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
  taxa.Name          = {}; % taxa name, lowest level.        Required field in WoRMS_AphiaID_taxa_match.m.
  taxa.Name_parent   = {}; % taxa parent name, if available. Required field in WoRMS_AphiaID_taxa_match.m.
  taxa.Name_original = {}; % name of original variable.      This field is not required for WoRMS_AphiaID_taxa_match.m.
  
  % Loop through and prase column names new fieldnames with associated ID
  for n_id = 1:numel(original_fields)
    taxa.ID{n_id}   = [num2str(n_id) 'id'];
    taxa.Name_original{n_id} = original_fields{n_id};
    % First, parse ecotaxa name to pull out taxa name and parent name
    names = strsplit(original_fields{n_id},'(');
    % Save name to table
    taxa.Name{n_id} = names{1};
    % See if parent name available
    if numel(names) > 1
      taxa.Name_parent{n_id} = erase(names{2},')');
    else
      taxa.Name_parent{n_id} = '';
    end
  end
  % Call script to go through and query AphiaID matches
  check_children_manual = false; % true = checks full classifcation of parent if child doesn't have a match and asks user. false = skips this.
  taxa = WoRMS_AphiaID_taxa_match(taxa,check_children_manual,'ecotaxa');
  fprintf('Saving taxa matches to %s\n',[pwd filesep save_taxa_filename])
  save(save_taxa_filename,'taxa');
else % Load data instead of rerunning everything
  fprintf('Loading..%s\n',save_taxa_filename)
  load(save_taxa_filename);
end

%% Match unknown LIVING taxa to their parent
% If parent is unknown, but still in the living categories, match to
% Eukaryota as suggested in OCB manual. 
idx_not_living = find(strcmp(taxa.Name,'not-living'));
for n_id = 1:idx_not_living-1
  % First, check if AphiaID has already been matched
  if strcmp(taxa.AphiaID{n_id},'NULL')
    % Parent has AphiaID, but no match found for child
    if ~strcmp(taxa.AphiaID_parent{n_id},'NULL')
      % match unknown taxonomic names to their parent AphiaID, when possible.
      taxa.AphiaID{n_id}          = taxa.AphiaID_parent{n_id};
      taxa.scientificNameID{n_id} = ['urn:lsid:marinespecies.org:taxname:' taxa.AphiaID{n_id}];
      taxa.Comment{n_id} = 'AphiaID --> AphiaID_parent';
      % Also set AphiaID_parent further along
      idx_parents = strcmp(taxa.Name_parent,taxa.Name{n_id});
      if sum(idx_parents) > 0
        taxa.AphiaID_parent(idx_parents) = taxa.AphiaID(n_id);
      end
      % Parent & child do not have a match but still living category
      % Match with Eukaryota
    elseif strcmp(taxa.AphiaID_parent{n_id},'NULL')
      taxa.scientificNameID{n_id} = 'urn:lsid:algaebase.org:taxname:86701';
      taxa.Comment{n_id} = 'miscellaneous living --> Eukaryota';
    end % IF parent has, or does not have, matched AphiaID
  end % IF AphiaID has already been matched
end

%% Non-conforming ROI's that need to be defined in a YAML-formatted namespace file
% These IDs (e.g. in our case taxa.Name{n_id}) are paired with definitions
% and are stored in a YAML  formatted file in order to serve as a
% machine-readable configuration file for anyone working with the data
% files.
for n_id = idx_not_living:size(taxa,1)
  % print to screen formatted text to insert into YAML file
  fprintf('  - id: %s\n',taxa.Name{n_id})
  fprintf('    definition: \n');
  % See if matches directly to OCB PTWG ids already defined
  if ismember(taxa.Name{n_id},ptwg)
    fprintf('    associated_terms: \n');
    fprintf('    - id: %s\n',['ptwg:' taxa.Name{n_id}]);
  elseif isfield(our_term,taxa.Name{n_id})
    fprintf('    associated_terms: \n');
    fprintf('    - id: %s\n',our_term.(taxa.Name{n_id}));
  else
    %fprintf('    associated_terms: \n');
    %fprintf('    - id: \n');
  end
  taxa.scientificNameID{n_id} = [namespace ':' taxa.Name{n_id}];
end
keyboard

%% Configure fieldnames to fit expected SeaBASS fields 
%Correct some variable names to match with SEABASS requirements
cols2 = strrep(cols2,'Site','station');
cols2 = strrep(cols2,'Station','station_alt_id');
cols2 = strrep(cols2,'Latitude_[degrees_north]','lat');
cols2 = strrep(cols2,'Longitude_[degrees_east]','lon');
% Variable SHOULD be bin_depth but fcheck throws and error
% "Measurement depth required as either a header value or a column in the
%  data block.  All data MUST have an associated depth"
%cols2= strrep(cols2,'Depth','bin_depth'); % bin_depth = Nominal or center depth for each data bin
cols2 = strrep(cols2,'Depth_[m]','depth');  % depth = Depth of measurement
cols2 = strrep(cols2,'Sampled_volume_[L]','volume');
% Set up units
units = cell(size(cols2));
if ismember('R2R_Event',cols2)
  units(1:9) = {'none' 'none' 'yyyymmdd' 'hh:mm:ss' 'none' 'degrees' 'degrees' 'm' 'L'};
  descr(1:9) = {'Station' 'Profile name' 'Date' 'Time (UTC)' 'Rolling Deck to Repository (R2R) Unique sample event number' 'Latitude in decimal degrees north' 'Longitude in decimal degrees north' 'Depth, midpoint of 5m depth bin' 'Total volume of water sampled within each depth bin'};
else
  units(1:8) = {'none' 'none' 'yyyymmdd' 'hh:mm:ss' 'degrees' 'degrees' 'm' 'L'};
  descr(1:8) = {'Station' 'Profile name' 'Date' 'Time (UTC)' 'Latitude in decimal degrees north' 'Longitude in decimal degrees north' 'Depth, midpoint of 5m depth bin' 'Total volume of water sampled within each depth bin'};
end

% Pull out abundance data [#/m^3]
idx_abun = contains(cols2,'_[#/m^3]');
cols2(idx_abun) = strcat('abun_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
units(idx_abun) = {'number/m^3'};
abun_desc = strcat('Abundance of',{' '},taxa.Name_original);
descr(idx_abun) = abun_desc; %strcat(abun_desc,{' '},erase(cols2(idx_abun),'abundance_'));

idx_biovol = contains(cols2,'_biovolume_[mm^3/L]');
cols2(idx_biovol) = strcat('biovol_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
% units(idx_biovol) = {'mm^3/L'};
% convert biovolume from [mm^3/L] to [µL/m^3] --> VSD * 1000.0; (1uL=1mm^3, 1000L=1m^3)
units(idx_biovol) = {'uL/m^3'};
biovol = table2array(odv2(:,idx_biovol)); % first need to convert to an array so can use the * operator
odv2(:,idx_biovol) = array2table(biovol.*1000.0); % [µL/m^3]
biovol_desc = strcat('Biovolume of',{' '},taxa.Name_original);
descr(idx_biovol) = biovol_desc; %strcat(biovol_desc,{' '},erase(cols2(idx_biovol),'biovolume_'));


% pull out esd data [mm]
idx_esd = contains(cols2,'_avgesd_[mm]');
cols2(idx_esd) = strcat('avgesd_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
units(idx_esd) = {'mm'};
esd_desc = strcat('Average equivalent spherical diameter of',{' '},taxa.Name_original);
descr(idx_esd) = esd_desc; %strcat(esd_desc,{' '},erase(cols2(idx_esd),'biovolume_'));

% Combine all units into a single string
colunits = strjoin(units,','); %[colstr,cols2{i},','];
% Combine column names into a single string
colstr = strjoin(cols2,','); % colstr=[colstr,cols2{i},','];

%% Generate table that includes all variables with description, units, size bin, etc
fields = table();
fields.VariableName = cols2';
fields.Units = units';
fields.Description = descr';
fields.Ecotaxa_original_name = repmat({''},size(fields.Units));
fields.Ecotaxa_original_name(idx_abun)   = taxa.Name_original;
fields.Ecotaxa_original_name(idx_biovol) = taxa.Name_original;
fields.Ecotaxa_original_name(idx_esd)    = taxa.Name_original;

% Add AphiaID
fields.WoRMS_AphiaID        = repmat({''},size(fields.Units));
fields.WoRMS_AphiaID_Parent = repmat({''},size(fields.Units));
fields.scientificNameID     = repmat({''},size(fields.Units));
fields.WoRMS_AphiaID(idx_abun)   = taxa.AphiaID;
fields.WoRMS_AphiaID(idx_biovol) = taxa.AphiaID;
fields.WoRMS_AphiaID(idx_esd)    = taxa.AphiaID;
fields.WoRMS_AphiaID_Parent(idx_abun)   = taxa.AphiaID_parent;
fields.WoRMS_AphiaID_Parent(idx_biovol) = taxa.AphiaID_parent;
fields.WoRMS_AphiaID_Parent(idx_esd)    = taxa.AphiaID_parent;
fields.scientificNameID(idx_abun)   = taxa.scientificNameID;
fields.scientificNameID(idx_biovol) = taxa.scientificNameID;
fields.scientificNameID(idx_esd)    = taxa.scientificNameID;
% remove NULL entries to avoid confusion
fields.WoRMS_AphiaID(strcmp(fields.WoRMS_AphiaID,'NULL')) = {''};
fields.WoRMS_AphiaID_Parent(strcmp(fields.WoRMS_AphiaID_Parent,'NULL')) = {''};

%% Write ParameterDescriptions table
%writetable(fields,fullfile('submit','UVP_ZOO_ParameterDescriptions.xlsx'))

%% Generate SeaBASS header text
% pull out max/min latitude, longitude, date, and time
latmax = max(odv2.latitude);
latmin = min(odv2.latitude);
lonmax = max(odv2.longitude);
lonmin = min(odv2.longitude);

% Find the first and last date/time for SEABASS metadata header
% since odv2 table has already been sorted by datetime (early on), just
% pull out first and last
datemax = odv2.date(end,:);
datemin = odv2.date(1  ,:);
timemax = odv2.time(end,:);
timemin = odv2.time(1  ,:);

% First - grab structure with basic metadata from SeaBASS_metadata
[~,original_file,ext] = fileparts(raw_rfile);
original_file = [original_file ext];
hdr = SeaBASS_metadata_headers(cruiseid,original_file);

% Write formatted header for sb file
hdr_SEABASS = format_UVP_SeaBASS_metadata_header;
                
%% Insert taxa definitions to comments
% AphiaID returns NULL if no match found, so replace all NaNs with NULL
% actually, use 'none' instead of 'NULL' so more human readable
taxa.AphiaID(strcmp(taxa.AphiaID,'NULL')) = {'none'};
taxa.AphiaID_parent(strcmp(taxa.AphiaID_parent,'NULL')) = {'none'};

% First, build text cells to insert
id_fields_definitions = '';
insert_ID_desc = cell(size(taxa,1)+2,1); % +2 because want to add note and column names
insert_ID_desc{1} = '! Taxomonic identifications have been matched to WoRMS AphiaIDs when possible.';
insert_ID_desc{2} = '!  ID    scientificName     AphiaID  data_provider_category          scientificNameID';
for nn = 1:size(taxa,1)
  insert_ID_desc{nn+2} = ['!  ' pad(taxa.ID{nn},5) ' ' pad(taxa.Name{nn},18) ' ' pad(taxa.AphiaID{nn},8) ' ' pad(taxa.Name_original{nn},31) ' ' taxa.scientificNameID{nn}];
  id_fields_definitions = [id_fields_definitions [taxa.ID{nn} ':' taxa.Name{nn}] ','];
end
% remove the last comma from the metadata header /id_fields_definitions
id_fields_definitions(end) = [];
id_fields_definitions = ['/id_fields_definitions=' id_fields_definitions];

idx_fields = find(contains(hdr_SEABASS,'/fields='));
hdr_SEABASS = [hdr_SEABASS(1:idx_fields-1); insert_ID_desc; '!'; id_fields_definitions; hdr_SEABASS(idx_fields:end)];

%% Specify file writing format
col_n = length(cols2);
col_n_text = find(contains(cols2,'lat')) - 1; %Latitude is the first DOUBLE column,
fmt_txt = repmat( '%s,', [1 col_n_text] );
fmt_dbl = repmat( '%f,', [1 (col_n - col_n_text)]);
% replace delimiter with new line at end of format string
fmt_dbl(end) = []; fmt_dbl = [fmt_dbl '\n'];
fmt = [fmt_txt fmt_dbl];
% reformat odv2 table to temporary variable to enable simply write to file
raw_write = table2cell(odv2); % convert to cell array to handle different variable types
raw_write = raw_write';            % transpose because fprintf function prints data columnwise
% convert NaNs to missing identifier
% raw_write(cellfun(@(x) isnumeric(x) && isnan(x), raw_write)) = {''}; % convert NaNs to '' for raw format
raw_write(cellfun(@(x) isnumeric(x) && isnan(x), raw_write)) = {-9999}; % missing=-9999 SeaBASS required header field

%% Write raw table to file
fprintf('\n  Deleting data file, if it exist (if it does not exist, you may get a WARNING, but this is OK: %s\n',raw_wfile);
eval(['delete ' raw_wfile])

fprintf('\n  Writing data in raw format to: %s\n',fullfile('submit',raw_wfile));
fileID = fopen(fullfile('submit',raw_wfile),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile('submit',raw_wfile))
  keyboard
end
fprintf(fileID,'%s\n',hdr_SEABASS{:});% write header
fprintf(fileID,fmt,raw_write{:});     % write data
fclose(fileID);                       % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
%fprintf('%s\n',hdr_SEABASS{:}) % write header
%fprintf(fmt,raw_write{:})      % write data

%% END OF MAIN FUNCTION - SUBFUNCTIONS BELOW
fprintf('workflow is finished\n')

%% FUNCTION META = FORMAT_UVP_SEABASS_METADATA_HEADER
  function hdr_SEABASS = format_UVP_SeaBASS_metadata_header
    % Needs hdr structure as defined in SeaBASS_metadata_headers
    if ~exist('hdr','var')
      error('Need "hdr" structure - called from SeaBASS_metadata_headers.m\n')
    end
    [~,original_file,ext] = fileparts(raw_rfile);
    original_file = [original_file ext];
    hdr_SEABASS={'/begin_header';
      ['/investigators='  hdr.investigators];
      ['/affiliations='   hdr.affiliations];
      ['/contact='        hdr.contact];
      ['/experiment='     hdr.experiment];
      ['/cruise='         hdr.cruise];
      ['/station='        hdr.station];
      ['/data_file_name=' raw_wfile]; 
      ['/original_file_name=' original_file]; % The original name of the data file, if different from the current /data_file_name. Designed to be a reference for the contributor.
      ['/documents='    strcat(strjoin(hdr.documents,','), [',' namespace '.yaml'])];
      ['/data_type='    hdr.data_type];
      ['/data_status='  hdr.data_status];
      ['/start_date='   datemin];
      ['/end_date='     datemax];
      ['/start_time='   timemin '[GMT]'];
      ['/end_time='     timemax '[GMT]'];
      ['/north_latitude=' num2str(latmax) '[DEG]'];
      ['/south_latitude=' num2str(latmin) '[DEG]'];
      ['/east_longitude=' num2str(lonmax) '[DEG]'];
      ['/west_longitude=' num2str(lonmin) '[DEG]'];
      '/water_depth=NA';
      ['/missing='    hdr.missing];
      ['/delimiter='  hdr.delimiter];
      ['/instrument_manufacturer='  hdr.inst_mfr];
      ['/instrument_model='         hdr.inst_model];
      ['/calibration_files='        hdr.calfiles];
      ['/calibration_date='         hdr.caldates];
      '!';};
    % Insert comments, then finish with /fields and /units
    hdr_SEABASS = [hdr_SEABASS; hdr.comments;...
                  '!';
                  ['/fields=' colstr];
                  ['/units=' colunits];
                  '/end_header'];
    % check if there is whitespace in any metadata headers
    % whitespace in comments is OKAY
    if any(contains(hdr_SEABASS,' ') & ~contains(hdr_SEABASS,'!'))
      fprintf('White space in metadata header, must remove to pass fcheck\n')
      keyboard
    end
  end %% FUNCTION META = FORMAT_UVP_SEABASS_METADATA_HEADER
end %% MAIN FUNCTION








