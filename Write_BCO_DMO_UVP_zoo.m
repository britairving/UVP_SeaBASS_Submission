function Write_BCO_DMO_UVP_zoo
%% FUNCTION WRITE_BCO_DMO_UVP_ZOO
% Description:
%   Reads in detailed ODV formatted file exported from Ecotaxa Particle
%   Module and creates a file for submission to BCO-DMO.
% 
% Steps:
%   1. Read in detailed ODV PAR file
%   2. Format fields and column names following BCO-DMO format
%         <https://www.bco-dmo.org/how-get-started>
%         BCO-DMO accepts all formats, so instead just change units for
%         clarity. E.g. [#_l-1] --> [#/L], [mm3_l-1] --> [mm^3/L]
%   3. Fill empty metadata cells with values from above
%   4. Delete size bins below UVP detection limit (~100 micrometers)
%
% References
%  https://www.bco-dmo.org/how-get-started#register_data
%  https://www.bco-dmo.org/submitting-data-spreadsheet
%
% Author:
%  Brita Irving <bkirving@alaska.edu> 
%  Andrew McDonnell <amcdonnell@alaska.edu>
%%
use_ids_for_fieldnames = 1; % instead of changing fields to abun_taxaname(parent_taxaname), just use abun_id#
%% Define read and write filenames
cruiseid   = 'RB1503';
projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);
odv_rfile          = fullfile(projectdir,'export_detailed_20210325_18_53','export_detailed_20210325_18_53_ZOO_odv.txt');
validated_profiles = fullfile(projectdir,'uvp5_sn009_2015_p16n_fully_validated_stations.txt');


%% Read in basic metadata for project
[~,original_file,ext] = fileparts(odv_rfile);
original_file = [original_file ext];
try
  % Change directory to project folder so can easily call metadata file
  pwd_now = pwd;
  cd(projectdir)
  eval(['hdr = ' cruiseid '_UVP_metadata(original_file)'])
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
odv_hdr = cell(20,1);
fileID = fopen(odv_rfile);
for nline = 1:20
  odv_hdr{nline} = fgetl(fileID);
end
fclose(fileID);
column_hdrline = find(contains(odv_hdr,'Cruise'));
% resize header
odv_hdr = odv_hdr(1:column_hdrline);
cols = strsplit(odv_hdr{column_hdrline},table_options.Delimiter); % split odv_hdr by odv delimiter
% store original column names
cols_orig = cols;
for ic = 1:numel(cols)
  cols{ic} = strrep(cols{ic},' ','_');
  cols{ic} = strrep(cols{ic},'__','_'); % clean up unnecessary underscores
end

% ecotaxa changed the format of the unit part of column names
% this happend sometime in 2019 (I think)
%       LPM       LPM_biovolume
% old: [#/m3]     [ppm]
% new: [# m-3]    [mm3 l-1]

% For simplicity, change column names back to original format
cols = strrep(cols,'_[#_m-3]','_[#/m^3]');
cols = strrep(cols,'_[mm3_l-1]','_[mm^3/L]');

%% Reformat variable names
meta_idx = find(contains(cols,'METAVAR'));

% remove unnecessary text from table column names
remove_text_from_columns = {'_METAVAR_TEXT_40','_METAVAR_TEXT_20', '_METAVAR_TEXT_10', '__METAVAR_DOUBLE','__PRIMARYVAR_DOUBLE','_degrees_east','_degrees_north'};
table_options.VariableNames = eras(table_options.VariableNames,remove_text_from_columns);
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

%% Remove unwanted columns
% define variables that we want to remove
remove_fields = {'DataOwner' 'Cruise' 'Instrument' 'SN' 'CTDrosette' 'extrames'};
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
  odv2 = odv2(:,new_order);
  % delete redundant variables that result from reordering
  odv2.date_1 = [];
  odv2.time_1 = [];
end
% sort by datetime then remove the datetime column
odv2.datetime = [];
cols2(i_datetime) = [];
cols2(end-1:end) = []; % date and time at the end


%% Query WoRMS to pull out AphiaID match for each taxonomic name
% save_taxa_filename = [cruiseid 'taxa_' date '.mat'];
save_taxa_filename = fullfile(projectdir,[cruiseid 'taxa.mat']);
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
  check_children_manual = false;% true = checks full classifcation of parent if child doesn't have a match and asks user. false = skips this.
  taxa = WoRMS_AphiaID_taxa_match(taxa,check_children_manual,'ecotaxa');
  fprintf('Saving taxa matches to %s\n',[pwd filesep save_taxa_filename])
  save(save_taxa_filename,'taxa');
else % Load data instead of rerunning everything
  fprintf('Loading... %s\n',save_taxa_filename)
  load(save_taxa_filename);
end

%% Configure fieldnames to align with SeaBASS fields 
%Correct some variable names to match with SEABASS requirements
cols2 = strrep(cols2,'Site','station');
cols2 = strrep(cols2,'Station','profile');
cols2 = strrep(cols2,'Latitude_[degrees_north]','lat');
cols2 = strrep(cols2,'Longitude_[degrees_east]','lon');
cols2 = strrep(cols2,'Rawfilename','UVP_rawfilename');
cols2 = strrep(cols2,'Depth_[m]','bin_depth');   
cols2 = strrep(cols2,'Sampled_volume_[L]','volume');
% Set up units
units(1:9) = {'none' 'none' 'none' 'yyyymmdd' 'hh:mm:ss' 'degrees_north' 'degrees_east' 'm' 'L'};
% Add parameter descriptions  
descr(1:9) = {'Site' 'Profile name' 'Raw file name from UVP' 'Date' 'Time (UTC)' 'Latitude in decimal degrees north' 'Longitude in decimal degrees north' 'Depth, midpoint of 5m depth bin' 'Total volume of water sampled within each depth bin'};

% Pull out abundance data [#/m^3]
idx_abun = contains(cols2,'[#/m^3]');
if use_ids_for_fieldnames
   cols2(idx_abun) = strcat('abun_',taxa.ID');
   abun_desc = strcat('Abundance of',{' '},taxa.Name_original);
else
  cols2(idx_abun) = erase(cols2(idx_abun),'_[#/m^3]');           % erase units from variable name
  cols2(idx_abun) = strcat('abun_',cols2(idx_abun));
  abun_desc = {'Abundance'};
end
units(idx_abun) = {'number/m^3'};
% add descriptions
descr(idx_abun) = abun_desc; %strcat(abun_desc,{' '},erase(cols2(idx_abun),'abundance_'));

% Pull out biovolume data [#/m^3]
idx_biovol = contains(cols2,'[mm^3/L]');
if use_ids_for_fieldnames
   cols2(idx_biovol) = strcat('biovol_',taxa.ID');
   biovol_desc = strcat('Biovolume of',{' '},taxa.Name_original);
else
  cols2(idx_biovol) = erase(cols2(idx_biovol),'_biovolume_[mm^3/L]');
  cols2(idx_biovol) = strcat('biovol_',cols2(idx_biovol));
  biovol_desc = {'Biovolume'};
end
units(idx_biovol) = {'mm^3/L'};
% add descriptions
descr(idx_biovol) = biovol_desc; %strcat(biovol_desc,{' '},erase(cols2(idx_biovol),'biovolume_'));

% pull out esd data [mm]
idx_esd = contains(cols2,'_avgesd_[mm]');
if use_ids_for_fieldnames
   cols2(idx_esd) = strcat('avgesd_',taxa.ID');
   esd_desc = strcat('Average equivalent spherical diameter of',{' '},taxa.Name_original);
else
  cols2(idx_esd) = erase(cols2(idx_esd),'_avgesd_[mm]');
  cols2(idx_esd) = strcat('avgesd_',cols2(idx_esd));  % Build variables as abun_1id,abun_2id... etc
  esd_desc = {'Average equivalent spherical diameter'};
end
units(idx_esd) = {'mm'};
descr(idx_esd) = esd_desc; %strcat(esd_desc,{' '},erase(cols2(idx_esd),'biovolume_'));

% Combine all units into a single string
colunits = strjoin(units,',');
% Combine column names into a single string
colstr   = strjoin(cols2,','); 

%% Only keep profiles that were fully validated
try
  % Read validated profiles
  fileID = fopen(validated_profiles,'r');
  profiles = fgetl(fileID);
  profiles = strtrim(strsplit(profiles,','));
  fclose(fileID);
  odv3 = odv2(ismember(odv2.profile,profiles),:);
  % fprintf(' START HERE: THROW OUT EMPTY DEPTH BINS!!\n')
  % abun = table2array(odv3(:,idx_abun));
  % empty_profiles = all(isnan(abun),2);
catch
  % Skip this part if there is not a list of validated stations
  odv3 = odv2;
end

%% Generate table that includes all variables with description, units, size bin, etc
fields = table();
fields.VariableName = cols2';
fields.Units = units';
fields.Description = descr';
if use_ids_for_fieldnames
  fields.Ecotaxa_original_name = repmat({''},size(fields.Units));
  fields.Ecotaxa_original_name(idx_abun)   = taxa.Name_original;
  fields.Ecotaxa_original_name(idx_biovol) = taxa.Name_original;
  fields.Ecotaxa_original_name(idx_esd)    = taxa.Name_original;
end
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

%% Add more description to temporary categories, if possible
% Definitions in [cruiseid]_UVP_metadata.m.
% for example...
% hdr.temporary_category.t002 = 'very long fibers or diatom chains or tentacles';
% hdr.temporary_category.t003 = 'portion of possible cnidarian';
% hdr.temporary_category.t004 = 'circular particle with one or two dark spots, unsure if living or non-living';
try
  if isfield(hdr,'temporary_category')
    temp_fields = fieldnames(hdr.temporary_category);
    % loop through temporary categories defined and add to description
    for nt = 1:numel(temp_fields)
      temp_name = temp_fields{nt};
      % Only replace the "temporary" text if there is text to replace it with
      if ~strcmp(hdr.temporary_category.(temp_name),'')
        replace_text = [temp_name '(temporary)'];
        idx_tempcat  = contains(fields.Description,replace_text);
        fields.Description(idx_tempcat) = strrep(fields.Description(idx_tempcat),'(temporary)',['(' hdr.temporary_category.(temp_name) ')']);
      end
    end
  end
catch
  fprintf('adding more descriptive text to the temporary categories did not work... \n')
end

%% Write ParameterDescriptions table
fprintf('Writing parameter descriptions to %s\n',fullfile(projectdir,'UVP_ZOO_ParameterDescriptions.xlsx'))
writetable(fields,fullfile(projectdir,'UVP_ZOO_ParameterDescriptions.xlsx'))

%% Generate header text to write to file
% first pull out max/min latitude, longitude, date, and time
latmax = max(odv3.latitude);
latmin = min(odv3.latitude);
lonmax = max(odv3.longitude);
lonmin = min(odv3.longitude);

% Find the first and last date/time for SEABASS metadata header
% since odv3 table has already been sorted by datetime (early on), just
% pull out first and last
datemax = odv3.date(end,:); 
datemin = odv3.date(1  ,:); 
timemax = odv3.time(end,:); 
timemin = odv3.time(1  ,:); 

% Write formatted header similar to SeaBASS format
hdr_text={'#begin_header';
  ['#investigators='  hdr.investigators];
  ['#affiliations='   hdr.affiliations];
  ['#contact='        hdr.contact];
  ['#experiment='     hdr.experiment];
  ['#cruise='         hdr.cruise];
  ['#data_type='      hdr.data_type];
  ['#start_date='     datemin];
  ['#end_date='       datemax];
  ['#start_time='     timemin '[GMT]'];
  ['#end_time='       timemax '[GMT]'];
  ['#north_latitude=' num2str(latmax) '[DEG]'];
  ['#south_latitude=' num2str(latmin) '[DEG]'];
  ['#east_longitude=' num2str(lonmax) '[DEG]'];
  ['#west_longitude=' num2str(lonmin) '[DEG]'];
  ['#missing_data='   hdr.missing];
  ['#instrument_manufacturer='  hdr.inst_mfr];
  ['#instrument_model='         hdr.inst_model];
  '#';};
% Insert comments
hdr_text = [hdr_text; hdr.comments.ZOO; '#end_header'];

%% Specify file writing format
if ~iscell(odv3.site)
  odv3.site = num2str(odv3.site);
end

col_n = length(cols2);
col_n_text = find(strcmp(cols2,'lat')) - 1; %Latitude is the first DOUBLE column
fmt_txt = repmat( '%s,', [1 col_n_text] );
fmt_dbl = repmat( '%f,', [1 (col_n - col_n_text)]);
% replace delimiter with new line at end of format string
fmt_dbl(end) = []; fmt_dbl = [fmt_dbl '\n'];
fmt = [fmt_txt fmt_dbl];
% reformat odv3 table to temporary variable to enable simply write to file

odv_write = table2cell(odv3); % convert to cell array to handle different variable types
odv_write = odv_write';      % transpose because fprintf function prints data columnwise
% convert NaNs to missing identifier
% https://www.bco-dmo.org/submitting-data-spreadsheet
% 3. Blanks vs. 0: There is a difference between a blank cell and a cell
%   with a value of zero (0). A blank cell denotes a missing or undefined
%   value, whereas a value of zero (0) means the value was measured as zero
%  (0). We usually change blank entries to the value of “nd” meaning “no
%  data”. (The quotes are not included.)

odv_write(cellfun(@(x) isnumeric(x) && isnan(x), odv_write)) = {-9999}; % missing=-9999 SeaBASS required header field

 %% Write odv table to file
% specify write filename
try
  wfile_name = fullfile(projectdir,[hdr.ecotaxa_name '_BCO-DMO_ZOO.txt']);
catch
  wfile_name = fullfile(projectdir,strrep(odv_rfile,'.txt','_BCO-DMO_ZOO.txt'));
end
% Delete file if it is release R0 - i.e. preliminary
if exist(wfile_name,'file') 
  delete(wfile_name)
end 

fprintf('\n  Writing data to: %s\n',wfile_name);
fileID = fopen(wfile_name,'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',wfile_name)
  keyboard
end
% write header
fprintf(fileID,'%s\n',hdr_text{:});
% Write column names
fprintf(fileID,'%s\n',colstr);
% Write column units
fprintf(fileID,'%s\n',colunits);
% write data
fprintf(fileID,fmt,odv_write{:});   
% close file
fclose(fileID);                     

fprintf('workflow is finished\n')

end %% MAIN FUNCTION WRITE_BCO_DMO_UVP_ZOO








