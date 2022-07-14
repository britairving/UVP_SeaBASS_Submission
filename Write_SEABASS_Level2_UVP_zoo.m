function Write_SEABASS_Level2_UVP_zoo
%% FUNCTION WRITE_SEABASS_LEVEL2_UVP_ZOO
% Description:
%   Reads in detailed ODV formatted file exported from Ecotaxa Particle
%   Module and creates a file for submission to NASA's SEABASS database
%   https://seabass.gsfc.nasa.gov/wiki/Data_Submission/
%
%   This will be Level 2 data (abundance and biovolumes)
%
% Steps:
%   1. USER writes metadata m file to define SeaBASS required headers
%      format: {cruiseid}_UVP_metadata.m
%   2. USER defines projectdir, cruiseid, odv_rfile (read file), sb_wfile 
%      (write file), and r2r_elog (if available and must edit 
%      UVP_merge_R2R.m for cruise specific r2r log formatting).
%   3. Script reads in detailed ODV ZOO file
%   4. Script formats fields and column names following SeaBASS format
%   5. Script queries WoRMS database for taxa scientificNames and scientificNameIDs
%   6. Script formats metadata header that will be printed to .sb file
%   7. Script writes sb file with all taxonomic data
%
% References:
%  Picheral M, Colin S, Irisson J-O (2017). EcoTaxa, a tool for the
%  taxonomic classification of images. http://ecotaxa.obs-vlfr.fr.
%
%  Neeley, A., S. Beaulieu, C. Proctor, I. Cetinic, J. Futrelle, I.
%  Ramos-Santos, H. Sosik, E. Devred, L. Karp-Boss, M. Picheral, N.
%  Poulton, C. Roesler, and A. Shepherd. (2021) Standards and practices for
%  reporting plankton and other particle observations from images. 38 pp.
%  DOI:10.1575/1912/27377
%
% Code availabile: https://github.com/britairving/UVP_submission_formatting
%
% Author:
%  Brita Irving <bkirving@alaska.edu>
%  with input from 
%     Andrew McDonnell <amcdonnell@alaska.edu>
%     Emmanuel Boss    <emmanuel.boss@maine.edu>
%     Lee Karp-Boss    <lee.karp-boss@maine.edu>
%% Uncertainty description
cfg.include_uncertainty_description = 1; % 1 = uncertainty description written to file header
cfg.remove_temporary_fields         = 0; % 1 = removes all temporary fields from odv data before writting

%% Define read and write filenames
%% EXPORTSNA cruise on the Sarmiento De Gamboa UVP6-LP deployed on a float
cruiseid  = 'SG2105_UVP6'; 
odv_rfile = fullfile('export_detailed_20220628_22_58','export_detailed_20220628_22_58_ZOO_odv.txt');
sb_wfile  = 'EXPORTS-EXPORTSNA_UVP6-TaxonomicLevel2_sdg_20210504-20210519_R1.sb';
r2r_elog  = 'EXPORTSNA_SarmientoDeGamboa_r2r_logs_original.xlsx';

%% EXPORTSNP survey cruise R/V Sally Ride
% cruiseid = 'SR1812';
% odv_rfile = fullfile('export_detailed_20210201_19_26','export_detailed_20210201_19_26_ZOO_odv.txt');
% sb_wfile = 'EXPORTS-EXPORTSNP_UVP5-TaxonomicLevel2_survey_20180814-20180909_R0.sb';

%% EXPORTSNP process cruise R/V Roger Revelle
% cruiseid = 'RR1813';
% odv_rfile = fullfile('export_detailed_20210128_08_11','export_detailed_20210128_08_11_ZOO_odv.txt');
% sb_wfile = 'EXPORTS-EXPORTSNP_UVP5-TaxonomicLevel2_process_20180814-20180909_R0.sb';


%% USER INPUT REQUIRED: Specify full path
if ismac 
  projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);
else
  projectdir = fullfile('D:','MATLAB','UVP_project_data',cruiseid);
end
odv_rfile = fullfile(projectdir,odv_rfile);
r2r_elog  = fullfile(projectdir,r2r_elog);


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
% Delete file if it is release R0 - i.e. preliminary
if exist(fullfile(submit_dir,sb_wfile),'file') && contains(sb_wfile,'R0')
  delete(fullfile(submit_dir,sb_wfile))
else % rename if file already exists so don't overwrite
  release_num = 1;
  while exist(fullfile(submit_dir,sb_wfile),'file')
    sb_wfile =  strrep(sb_wfile,['R' num2str(release_num) '.sb'],['R' num2str(release_num+1) '.sb']);
    release_num = release_num + 1;
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

%% Read ODV zoo data
fprintf('\n  Reading ODV Zooplankton data from... %s\n',odv_rfile)
% check that file exists
if ~exist(odv_rfile,'file')
  fprintf('File does not exist - make sure path and name are correct: %s\n',odv_rfile)
end
% STEP1 - read header and table options
% Get file options delimiter can be ';' or '\t' (tab)
table_options = detectImportOptions(odv_rfile);
% Read header
% column names usually at line 7, but go to line 20 in case extra
% header information is in file
hdr_odv = cell(20,1);
fileID = fopen(odv_rfile);
for nline = 1:20
  hdr_odv{nline} = fgetl(fileID);
end
fclose(fileID);
column_hdrline = find(contains(hdr_odv,'Cruise'));
% resize header
hdr_odv = hdr_odv(1:column_hdrline);
cols = strsplit(hdr_odv{column_hdrline},table_options.Delimiter); % split hdr_odv by odv delimiter
% store original column names
cols_orig = cols;
for ic = 1:numel(cols)
  cols{ic} = strrep(cols{ic},' ','_');
  cols{ic} = strrep(cols{ic},'__','_'); % clean up unnecessary underscores
end

% ecotaxa changed the format of the unit part of column names
% this happend sometime in 2019 (I think)
% old: [#/m3]     [ppm]
% new: [# m-3]    [mm3 l-1]

% For simplicity, change column names back to original format
cols = strrep(cols,'_[#_m-3]','_[#/m^3]');
cols = strrep(cols,'_[mm3_l-1]','_[mm^3/L]');

%% Reformat variable names
meta_idx = find(contains(cols,'METAVAR'));

% remove unnecessary text from table column names
remove_text_from_columns = {'_METAVAR_TEXT_40','_METAVAR_TEXT_20', '_METAVAR_TEXT_10', '__METAVAR_DOUBLE','__PRIMARYVAR_DOUBLE','_degrees_east','_degrees_north'};
table_options.VariableNames = erase(table_options.VariableNames,remove_text_from_columns);
% now remove unnecessary text from original column names
remove_text_from_columns = strrep(remove_text_from_columns,'_',':');
remove_text_from_columns = strrep(remove_text_from_columns,'::',':');
cols = erase(cols,remove_text_from_columns);

% make sure all variables read in correctly as number and not strings
table_options.VariableTypes(meta_idx(end):end) = {'double'};
table_options.VariableNames(1:meta_idx(end))   = lower(table_options.VariableNames(1:meta_idx(end)));
table_options.VariableNames{contains(table_options.VariableNames,'yyyy_mm_ddhh_mm')} = 'datetime';
table_options.VariableNames{contains(table_options.VariableNames,'station')} = 'profile';
table_options.VariableNames{contains(table_options.VariableNames,'Depth_m')} = 'Depth';

% STEP2 - read ODV particle data from file
odv = readtable(odv_rfile,table_options);

%% Read R2R eventlog
try
  % Change directory to project folder so can easily call metadata file
  pwd_now = pwd;
  cd(projectdir)
  eval(['[odv,cols] = ' cruiseid '_UVP_R2R_merge(odv,cols);'])
  % Change directory back to previous working directory
  cd(pwd_now);
catch % catch and explain why script stopped
  fprintf('Cannot load R2R merge script for this project\n')
  %fprintf('Need to set up %s_UVP_R2R_merge.m script, see example provided\n',cruiseid)
  %odv.R2R_Event = [];
  %keyboard
end

%% Define new table
% this is just to preserve original data so any manipulations can be
% tracked
odv2  = odv;
cols2 = cols;

%% Fill empty meta variables
% This is necessary because ODV format has multiple rows of data for each
% profile but after the first row, all subsequent rows for that profile are
% empty
fnames = fieldnames(odv2);
for iv = 1:numel(meta_idx)
  vname = fnames{meta_idx(iv)};
  odv2.(vname) = fillmissing(odv2.(vname),'previous');
end % loop through meta variables

% also fill out r2r event
if ismember('R2R_Event',fnames)
  % first, need to replace 0x0 double entries with ''
  odv2.R2R_Event(cellfun('isempty', odv2.R2R_Event)) = {''};
  odv2.R2R_Event = fillmissing(odv2.R2R_Event,'previous');
end

%   %% Identify all empty variables - i.e. no real values
%   rm_fields = [];
%   for irm = meta_idx(end)+1:size(odv2,2)
%     nam = fnames{irm};
%     if all(isnan(odv2.(nam)))
%       rm_fields = [rm_fields; irm];
%     elseif all(odv2.(nam) == 0)
%       rm_fields = [rm_fields; irm];
%     end
%   end
%
%   %% Remove empty variables
%   odv2(:,rm_fields) = [];
%   cols2(rm_fields)  = [];

%% Remove temporary fields
if cfg.remove_temporary_fields
  temporary_fields = find(contains(cols2,'temporary'));
  odv2(:,temporary_fields) = [];
  cols2(temporary_fields)  = [];
end

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
  % delete redundant variables that result from reordering
  odv2.date_1 = [];
  odv2.time_1 = [];
end
% sort by datetime then remove the datetime column
odv2.datetime     = [];
cols2(i_datetime) = [];

cols2(end-1:end) = []; % date and time at the end

% Now move r2r event to just after depth
if ismember('R2R_Event',cols2)
  ncols = numel(odv2.Properties.VariableNames);
  % find index location of variables we want to reorder
  i_r2revent = find(strcmp(odv2.Properties.VariableNames,'R2R_Event'));
  i_time     = find(strcmp(odv2.Properties.VariableNames,'time'));
  new_order  = [1:i_time i_r2revent i_time+1:ncols];
  % Reorder cell array with full column names
  cols2 = cols2(new_order); % replaced cols2=[cols2(1,1:dt-1) , cols2(1,d) , cols2(1,d+1) , cols2(1,dt+1:d-1) , cols2(1,d+2:end)];
  odv2 = odv2(:,new_order);
  % delete redundant variables that result from reordering
  odv2.R2R_Event_1    = [];
  cols2(i_r2revent+1) = [];
end


%% Query WoRMS to pull out AphiaID match for each taxonomic name
% save_taxa_filename = [cruiseid 'taxa_' date '.mat'];
save_taxa_filename = fullfile(projectdir,[cruiseid 'taxa_level2.mat']);
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
  taxa.ID            = {}; % suffix added to abun_zoop_, biovol_, avgesd_ fields. This field is not required for WoRMS_AphiaID_taxa_match.m
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

%% Limit data to specific taxa as defined by cfg.limit_to_taxa
% Cell array with taxonomic names selected by the user
% For example:  cfg.limit_to_taxa = {'Salpida' 't008'};   
% Used for EXPORTSNP UVP Dataset:  salps and salp pellets is fine with me
% as it was carefully checked by experts (Deb, Karen).
if isfield(cfg,'limit_to_taxa')
  idx_taxa = ismember(taxa.Name,cfg.limit_to_taxa);
  taxa = taxa(idx_taxa,:);
  % Pull out abundance data [#/m^3],biovolume data [mm^3/L], & esd data [mm]
  idx_taxa = find(contains(cols2,{'_[#/m^3]' '_biovolume_[mm^3/L]' '_avgesd_[mm]'}));
  idx_remove_taxa = ~contains(cols2(idx_taxa),cfg.limit_to_taxa); % do not include metadata columns
  idx_remove_taxa = idx_taxa(idx_remove_taxa);            % reindex to full size
  cols2(idx_remove_taxa)  = [];
  odv2(:,idx_remove_taxa) = [];
  
  % Redo taxa ID, in case any were removed
  % First step | parse taxonomic names from Ecotaxa fields
  % Only need to look at abundance names, because column names are the same
  % for biovolume and avgesd.
  original_fields = cols2(contains(cols2,'[#/m^3]'));
  original_fields = erase(original_fields,'_[#/m^3]');
  for n_id = 1:numel(original_fields)
    taxa.ID{n_id}   = [num2str(n_id) 'id'];
  end
end

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
% Notes about variable names from NASA
%  I did a few modifications to the field names in the taxonomic file to
%  adjust to our system who can only handle one type of unit per field.
%  Abun was taken with cell/L units, so we created a new field call
%  abun_zoop, same for biovol.
% 
% For the avgsed, we went with a long torturing one! The IFCB had set
% equivalent_spherical_diameter for the level-1b in um units, so we created
% yours as equivalent_spherical_diameter_avgmm.

% Pull out abundance data [#/m^3]
idx_abun = contains(cols2,'_[#/m^3]');
cols2(idx_abun) = strcat('abun_zoop_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
units(idx_abun) = {'organisms/L'}; % NASA prefers 'organisms/L' to 'number/m^3'
abun = table2array(odv2(:,idx_abun)); % first need to convert to an array so can use the * operator
odv2(:,idx_abun) = array2table(abun./1000.0); % convert to organisms/L = number/m^3 / 1000.00.
abun_desc = strcat('Abundance of',{' '},taxa.Name_original);
descr(idx_abun) = abun_desc; %strcat(abun_desc,{' '},erase(cols2(idx_abun),'abundance_'));


idx_biovol = contains(cols2,'_biovolume_[mm^3/L]');
cols2(idx_biovol) = strcat('biovol_zoop_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
% units(idx_biovol) = {'mm^3/L'};
% convert biovolume from [mm^3/L] to [�L/m^3] --> VSD * 1000.0; (1uL=1mm^3, 1000L=1m^3)
units(idx_biovol) = {'uL/m^3'};
biovol = table2array(odv2(:,idx_biovol)); % first need to convert to an array so can use the * operator
odv2(:,idx_biovol) = array2table(biovol.*1000.0); % [�L/m^3]
biovol_desc = strcat('Biovolume of',{' '},taxa.Name_original);
descr(idx_biovol) = biovol_desc; %strcat(biovol_desc,{' '},erase(cols2(idx_biovol),'biovolume_'));


% pull out esd data [mm]
idx_esd = contains(cols2,'_avgesd_[mm]');
cols2(idx_esd) = strcat('equivalent_spherical_diameter_avgmm_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
units(idx_esd) = {'mm'};
esd_desc = strcat('Average equivalent spherical diameter of',{' '},taxa.Name_original);
descr(idx_esd) = esd_desc; %strcat(esd_desc,{' '},erase(cols2(idx_esd),'biovolume_'));


% Combine all units into a single string
colunits = strjoin(units,','); %[colstr,cols2{i},','];
% Combine column names into a single string
colstr = strjoin(cols2,','); % colstr=[colstr,cols2{i},','];

%% Generate SeaBASS header text
% Pull out raw filename
[~,original_file,ext] = fileparts(odv_rfile);
original_file = [original_file ext];

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

hdr_SEABASS={'/begin_header';
  ['/investigators='  hdr.investigators];
  ['/affiliations='   hdr.affiliations];
  ['/contact='        hdr.contact];
  ['/experiment='     hdr.experiment];
  ['/cruise='         hdr.cruise];
  ['/station='        hdr.station];
  ['/data_file_name=' sb_wfile];
  ['/associated_files=' original_file]; % The original name of the data file exported from Ecotaxa
  '/associated_file_types=raw';
  ['/documents='      strjoin(hdr.documents.ZOO,',')];
  ['/data_type='      hdr.data_type];
  ['/data_status='    hdr.data_status.ZOO];
  ['/start_date='     datemin];
  ['/end_date='       datemax];
  ['/start_time='     timemin '[GMT]'];
  ['/end_time='       timemax '[GMT]'];
  ['/north_latitude=' num2str(latmax) '[DEG]'];
  ['/south_latitude=' num2str(latmin) '[DEG]'];
  ['/east_longitude=' num2str(lonmax) '[DEG]'];
  ['/west_longitude=' num2str(lonmin) '[DEG]'];
  '/water_depth=NA';
  ['/missing='        hdr.missing];
  ['/delimiter='      hdr.delimiter];
  ['/instrument_manufacturer='  hdr.inst_mfr];
  ['/instrument_model='         hdr.inst_model];
  ['/calibration_files='        hdr.calfiles];
  ['/calibration_date='         hdr.caldates];
  '!'};
  if cfg.include_uncertainty_description
  hdr.comments.ZOO = [hdr.comments.ZOO; ...
    '!';...
    '! Uncertainties can be calculated based on counting statistics';...
    '! For example:' ;...
    '!  relative uncertainty [none]';...
    '!  relative_unc = sqrt(abun_zoop_{#}id*volume)/(abun_zoop_{#}id*volume)';
    '!  uncertainty in the abundance [organisms/L]';...
    '!  abun_zoop_{#}id_unc = relative_unc*abun_zoop_{#}id';
    '!  uncertainty in the biovolume [uL/m^3]';...
    '!  biovol_zoop_{#}id_unc = relative_unc*biovol_zoop_{#}id'];%...
  %'!  uncertainty in the average equivalent spherical diameter [mm]';...
  %'!  avgesd_id{#}_unc = ... not sure'};];
  end
  % Insert comments, then finish with /fields and /units
  hdr_SEABASS = [hdr_SEABASS; hdr.comments.ZOO;...
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
odv_write = table2cell(odv2); % convert to cell array to handle different variable types
odv_write = odv_write';            % transpose because fprintf function prints data columnwise
% convert NaNs to missing identifier
% odv_write(cellfun(@(x) isnumeric(x) && isnan(x), odv_write)) = {''}; % convert NaNs to '' for ODV format
odv_write(cellfun(@(x) isnumeric(x) && isnan(x), odv_write)) = {-9999}; % missing=-9999 SeaBASS required header field

%% Write odv table to file
fprintf('\n  Deleting data file, if it exist (if it does not exist, you may get a WARNING, but this is OK: %s\n',sb_wfile);
eval(['delete ' sb_wfile])

fprintf('\n  Writing data in ODV format to: %s\n',fullfile(submit_dir,sb_wfile));
fileID = fopen(fullfile(submit_dir,sb_wfile),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile(submit_dir,sb_wfile))
  keyboard
end
fprintf(fileID,'%s\n',hdr_SEABASS{:});% write header
fprintf(fileID,fmt,odv_write{:});     % write data
fclose(fileID);                       % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
%fprintf('%s\n',hdr_SEABASS{:}) % write header
%fprintf(fmt,odv_write{:})      % write data

%% END OF MAIN FUNCTION - SUBFUNCTIONS BELOW
fprintf('workflow is finished\n')

%% FUNCTION META = FORMAT_UVP_SEABASS_METADATA_HEADER
  function hdr_SEABASS = format_UVP_SeaBASS_metadata_header

  end %% FUNCTION META = FORMAT_UVP_SEABASS_METADATA_HEADER
end %% MAIN FUNCTION








