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
%   1. Read in detailed ODV ZOO file
%   2. Format fields and column names following SeaBASS format
%   3. Define metadata header that will be printed to .sb file
%   4. Write sb file with all taxonomic data
%
% References:
%  Picheral M, Colin S, Irisson J-O (2017). EcoTaxa, a tool for the
%  taxonomic classification of images. http://ecotaxa.obs-vlfr.fr.
%
% Author:
%  Brita Irving <bkirving@alaska.edu>
%  Andrew McDonnell <amcdonnell@alaska.edu>
%% Uncertainty description
include_uncertainty_description   = 1;
remove_all_but_salp_and_salpfeces = 1;
%% Define read and write filenames
%% EXPORTSNP survey cruise R/V Sally Ride
cruiseid = 'SR1812';
odv_rfile = 'D:\MATLAB\EXPORTS_UVP_SEABASS\export_detailed_20210219_19_10\export_detailed_20210219_19_10_ZOO_odv.txt';
odv_wfile = 'EXPORTS-EXPORTSNP_UVP5-Taxonomic_survey_20180814-20180909_R0.sb';
r2r_elog  = 'D:\MATLAB\EXPORTS_UVP_SEABASS\R2R_ELOG_SR1812_FINAL_EVENTLOG_20180913_022931.xlsx';
%% EXPORTSNP process cruise R/V Roger Revelle
% cruiseid = 'RR1813';
% odv_rfile = 'D:\MATLAB\EXPORTS_UVP_SEABASS\export_detailed_20210128_08_11\export_detailed_20210128_08_11_ZOO_odv.txt';
% odv_wfile = 'EXPORTS-EXPORTSNP_UVP5-Taxonomic_process_20180814-20180909_R0.sb';
% r2r_elog  = 'D:\MATLAB\EXPORTS_UVP_SEABASS\R2R_ELOG_RR1813_FINAL_EVENTLOG_20180912_170812.xlsx';

if remove_all_but_salp_and_salpfeces
  odv_wfile = strrep(odv_wfile,'Taxonomic','TaxonomicLevel2_Salps');
end
% Define namespace YAML-formatted filename used to define non-conforming
% Use the same namespace for both EXPORTSNP cruises RR1813 and
% SR1812 because submitted at the same time... February 2021
switch cruiseid
%   case {'SR1812' 'RR1813'}
%     namespace = 'EXPORTSNP_Ecotaxa';
  otherwise
    namespace = ['EXPORTS_' cruiseid];
end

% specify write file, this is not computer specific because uses same path
% Write files to submit folder
if ~isdir('submit')
  mkdir('submit');
end

% Delete file if it is release R0 - i.e. preliminary
if exist(fullfile('submit',odv_wfile),'file') && contains(odv_wfile,'R0')
  delete(fullfile('submit',odv_wfile))
else % rename if file already exists so don't overwrite
  release_num = 1;
  while exist(fullfile('submit',odv_wfile),'file')
    odv_wfile =  strrep(odv_wfile,['R' num2str(release_num) '.sb'],['R' num2str(release_num+1) '.sb']);
    release_num = release_num + 1;
  end
end

% Store non-conforming terms that the OCB PTWG (Ocean Carbon and
% Biogeochemistry Phytoplankton Taxonomy Working Group) created. Created to
% define several standardized names for common terms that are not currently
% defined by WoRMS or Algae Base.
% The term ‘other’ should only be used to describe a non-living particle.
ptwg = {'bad_image' 'bead' 'bubble' 'detritus' 'fecal_pellet' 'other'};
our_term.badfocus = 'ptwg:bad_image';
our_term.feces    = 'ptwg:fecal_pellet';

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
hdr = cell(20,1);
fileID = fopen(odv_rfile);
for nline = 1:20
  hdr{nline} = fgetl(fileID);
end
fclose(fileID);
column_hdrline = find(contains(hdr,'Cruise'));
% resize header
hdr = hdr(1:column_hdrline);
cols = strsplit(hdr{column_hdrline},table_options.Delimiter); % split hdr by odv delimiter
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
  r2r = readtable(r2r_elog,'FileType','spreadsheet');
  % Initialize r2r event field
  odv.R2R_Event = cell(size(odv.cruise));
  %odv.r2r_cast  = cell(size(odv.cruise));
  cols = [cols, 'R2R_Event'];
  % Pull out r2r event based on matching event
  switch cruiseid
    case 'SR1812'
      [rawfilenames,iu] = unique(odv.rawfilename);
      % remove empty entry (necessary because odv format only has single
      % entry per profile
      if isempty(rawfilenames{1})
        rawfilenames(1) = [];
        iu(1) = [];
      end
      % pull out profiles from uvp data
      profiles = erase(odv.profile(iu),{'ctd00' 'ctd0' 'ctd'});
      % loop through profiles/rawfilenames and pull out r2revent
      for nr = 1:numel(rawfilenames)
        % math by UVP rawfilename, then station or cast
        mtch = contains(r2r.Comment,rawfilenames{nr}) & ...
          ( strcmp(strrep(r2r.Station,' ','_'),odv.site{iu(nr)}) | strcmp(r2r.Cast,profiles{nr}) );
        if sum(mtch) == 1
          odv.R2R_Event{iu(nr)} = r2r.R2R_Event{mtch};
          %odv.r2r_cast{iu(nr)} = r2r.Cast{mtch};
        else
          fprintf('r2r event not found, stopping here\n')
          keyboard
        end
      end
    case 'RR1813'
      [sites,iu] = unique(odv.site);
      % remove empty entry (necessary because odv format only has single
      % entry per profile
      if isempty(sites{1})
        sites(1) = [];
        iu(1) = [];
      end
      % pull out profiles from uvp data
      profiles1 = erase(odv.profile(iu),{'export00' 'export0' 'export'});
      profiles2 = strcat('SIO_',pad(profiles1,3,'left','0'));
      % loop through profiles/rawfilenames and pull out r2revent
      for nr = 1:numel(sites)
        sdate = sites{nr}(1:8);
        % math by UVP rawfilename, then station or cast
        mtch = contains(r2r.Event,sites{nr}) | contains(r2r.R2R_Event,sites{nr});
        if sum(mtch) ~= 1
          mtch = contains(r2r.Event,sdate) & contains(r2r.Instrument,'CTD911') & ( strcmp(r2r.Cast,profiles1{nr}) | strcmp(r2r.Cast,profiles2{nr}) );
        end
        if sum(mtch) ~= 1
          % for most of the casts, r2r event log has casts as some
          % combindation of SIO and sequential cast number.
          mtch = contains(r2r.Event,sdate) & contains(r2r.Instrument,'CTD911') & contains(r2r.Cast,'SIO','IgnoreCase',true) & contains(r2r.Cast,profiles1{nr}) ;
        end
        % set r2r event
        odv.R2R_Event{iu(nr)} = r2r.R2R_Event{mtch};
        %odv.r2r_cast{iu(nr)}  = r2r.Cast{mtch};
      end
    otherwise
      fprintf('R2R event log not match up yet, ignoring for now\n')
      odv.R2R_Event = [];
  end
catch
  fprintf('Could not read R2R event log file, ignoring for now\n')
end

%% Define new table
% this is just to preserve original data so any manipulations can be
% tracked
odv2  = odv;
cols2 = cols;


%% Remove all except salp and salp pellets
if remove_all_but_salp_and_salpfeces
  % but need to make sure that we give similar IDs (e.g., for salp feces)
  % salps and salp pellets is fine with me as it was carefully checked by experts (Deb, Karen).
  % Perhaps in the meantime we submit a level 2 file only for salps and salp pellets?  That could be done rather quickly by using Brita's existing level 2 file and removing all non-salp fields.
  keep_taxa = {'Salpida' 't008'};
 
  % Pull out abundance data [#/m^3],biovolume data [mm^3/L], & esd data [mm]
  idx_taxa = find(contains(cols2,{'_[#/m^3]' '_biovolume_[mm^3/L]' '_avgesd_[mm]'}));
  
  idx_remove_taxa = ~contains(cols2(idx_taxa),keep_taxa); % do not include metadata columns
  idx_remove_taxa = idx_taxa(idx_remove_taxa);            % reindex to full size
  cols2(idx_remove_taxa)  = [];
  odv2(:,idx_remove_taxa) = [];
end

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
% % Decided NOT to remove temporary fields based on email correspondence
% with PIs.
%temporary_fields = find(contains(cols2,'temporary'));
%odv2(:,temporary_fields) = [];
%cols2(temporary_fields)  = [];

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

%% Non-conforming ROI's that need to be defined in a YAML-formatted namespace file
% These IDs (e.g. in our case taxa.Name{n_id}) are paired with definitions
% and are stored in a YAML  formatted file in order to serve as a
% machine-readable configuration file for anyone working with the data
% files.
idx_not_living = find(strcmp(taxa.Name,'not-living'));
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
% writetable(fields,fullfile('submit',[cruiseid '_UVP_Taxonomic_ParameterDescriptions.xlsx']))

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
[~,original_file,ext] = fileparts(odv_rfile);
original_file = [original_file ext];
hdr = SeaBASS_metadata_headers(cruiseid,original_file);

% Write formatted header for sb file
hdr_SEABASS = format_UVP_SeaBASS_metadata_header;
if remove_all_but_salp_and_salpfeces
  % Sally Ride
  expert_image_val = '!  Expert image validation: Debbie Steinberg and Karen S. Stamieszkin'; 
  idx_find = find(contains(hdr_SEABASS,'Image validation'));
  hdr_SEABASS = [hdr_SEABASS(1:idx_find); expert_image_val; hdr_SEABASS(idx_find+1:end)];
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
fprintf('\n  Deleting data file, if it exist (if it does not exist, you may get a WARNING, but this is OK: %s\n',odv_wfile);
eval(['delete ' odv_wfile])

fprintf('\n  Writing data in ODV format to: %s\n',fullfile('submit',odv_wfile));
fileID = fopen(fullfile('submit',odv_wfile),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile('submit',odv_wfile))
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
    % Needs hdr structure as defined in SeaBASS_metadata_headers
    if ~exist('hdr','var')
      error('Need "hdr" structure - called from SeaBASS_metadata_headers.m\n')
    end
    [~,original_file,ext] = fileparts(odv_rfile);
    original_file = [original_file ext];
    hdr_SEABASS={'/begin_header';
      ['/investigators='  hdr.investigators];
      ['/affiliations='   hdr.affiliations];
      ['/contact='        hdr.contact];
      ['/experiment='     hdr.experiment];
      ['/cruise='         hdr.cruise];
      ['/station='        hdr.station];
      ['/data_file_name=' odv_wfile]; 
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
    if include_uncertainty_description
      hdr.comments = [hdr.comments; ...
        '!';...
        '! Uncertainties can be calculated based on counting statistics';...
        '! For example:' ;...
        '!  relative uncertainty [none]';...
        '!  relative_unc = sqrt(abun_id{#}*volume)/(abun_id{#}*volume)';
        '!  uncertainty in the abundance [number/m^3]';...
        '!  abun_id{#}_unc = relative_unc*abun_id{#}';
        '!  uncertainty in the biovolume [mm^3/L]';...
        '!  biovol_id{#}_unc = relative_unc*biovol_id{#}'];%...
        %'!  uncertainty in the average equivalent spherical diameter [mm]';...
        %'!  avgesd_id{#}_unc = ... not sure'};];
    end
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








