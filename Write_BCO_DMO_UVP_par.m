function Write_BCO_DMO_UVP_par
%% FUNCTION WRITE_BCO_DMO_UVP_PAR
% Description:
%   Reads in detailed ODV formatted file exported from Ecotaxa Particle
%   Module and creates a file for submission to NASA's SEABASS database
%   https://seabass.gsfc.nasa.gov/wiki/Data_Submission/
% 
% Steps:
%   1. Read in detailed ODV PAR file
%   2. Format fields and column names following BCO-DMO format
%         <https://www.bco-dmo.org/how-get-started>
%         BCO-DMO accepts all formats, so instead just change units for
%         clarity. E.g. [#_l-1] --> [#/L], [mm3_l-1] --> [mm^3/L]
%   3. Fill empty metadata cells with values from above
%   4. Delete size bins below UVP detection limit (~100 micrometers)
%   5. 
%
% References
%  https://www.bco-dmo.org/how-get-started#register_data
%  https://www.bco-dmo.org/submitting-data-spreadsheet
%
% Author:
%  Brita Irving <bkirving@alaska.edu> 
%  Andrew McDonnell <amcdonnell@alaska.edu>
%% 0 | Set script flags
% See section 6 | DO NOT DELETE higher size bins after first size bin with real data
% delete_all_empty_size_bins == 1 | Deletes ALL empty size bins
% delete_all_empty_size_bins == 0 | Keeps all size bins after the first
% one that contains data, even if that size bin contains no real data. E.g.
% could be useful for provenance or data archival purposes. All size bins
% below the first one to contain data are deleted.
delete_all_empty_size_bins = 0; 
remove_ctd_fields = 1; 

%% Define read and write filenames
cruiseid   = 'RB1503';
projectdir = fullfile('/Users/bkirving/Documents/MATLAB/UVP_project_data',cruiseid);
odv_rfile  = fullfile(projectdir,'export_detailed_20210325_18_53','export_detailed_20210325_18_53_PAR_odv.txt');

 % specify write file, this is not computer specific because uses same path
% Write files to submit folder
submit_dir = fullfile(projectdir,'submit');
if ~isdir(submit_dir)
  mkdir(submit_dir);
end

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
%% Read ODV particle data
fprintf('\n  Reading ODV Particle data from %s\n',odv_rfile)
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
for ic = 1:numel(cols)
  cols{ic} = strrep(cols{ic},' ','_');
end

% ecotaxa changed the format of the unit part of column names
% this happend sometime in 2019 (I think)
% For simplicity, change column names back to original format

if any( contains(cols,'_[#_l-1]') )
  cols = strrep(cols,'_[#_l-1]','[#/L]');
end
if any (contains(cols,')_[mm3_l-1]') )
  cols = strrep(cols,')_[mm3_l-1]',')[mm^3/L]');
end
% Reformat variable names
meta_idx = find(contains(cols,'METAVAR'));

% remove unnecessary text from column names
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

%% Identify all empty variables for removal - i.e. no real values
rm_fields = [];
for irm = meta_idx(end)+1:size(odv2,2)
  nam = fnames{irm};
  if ~iscell(odv2.(nam))
    if all(isnan(odv2.(nam)))
      rm_fields = [rm_fields; irm];
    elseif all(odv2.(nam) == 0)
      rm_fields = [rm_fields; irm];
    end
  end
end

%% DO NOT DELETE higher size bins after first size bin with real data
if delete_all_empty_size_bins == 0
  % For example, in many UVP5 HD datasets, the 50.8-64µm contains data, but
  % the next size bin does not. For the non-HD units, the 102-128µm size bin
  % contains data but the 128-161µm does not.
  % Detection limit for UVP is approximately 100 µm but there is commonly
  % data in the 50.8-64µm size bin and we want to maintain all columns with
  % data for posterity sake.
  idx_nperL  = find(contains(cols,'[#/L]')); % Indices of number concentration variables
  idx_biovol = find(contains(cols,'[mm^3/L]')); % Indices of biovolume variables
  % Find first size class that contains data
  nperL_data     = table2array(odv(:,idx_nperL));
  idx_firstgood  = find(any(isfinite(nperL_data),1),1);
  idx_keep_nperL = idx_nperL(idx_firstgood):idx_nperL(end);
  % Find first size class that contains data
  biovol_data    = table2array(odv(:,idx_biovol));
  % Biovolume fields with no value are exported from ecotaxa as 0, instead of
  % blank. Therefore, replace all biovolume == 0 with NaN to avoid confusion.
  biovol_data(biovol_data == 0) = NaN;
  idx_firstgood   = find(any(isfinite(biovol_data),1),1);
  idx_keep_biovol = idx_biovol(idx_firstgood):idx_biovol(end);
  ignore_empty    = ismember(rm_fields, [idx_keep_nperL idx_keep_biovol] );
  rm_fields(ignore_empty) = []; % keep the columns meeting the above criteria
end

%% Identify all fields after last biovolume field for removal - i.e. any CTD fields or
% ancillary variables
if remove_ctd_fields
  rm_fields_after_biovol = [idx_keep_biovol(end)+1:numel(cols2)]';
  ignore_empty = ismember(rm_fields_after_biovol,rm_fields);
  rm_fields_after_biovol(ignore_empty) = [];
  rm_fields = [rm_fields; rm_fields_after_biovol];
end

%% Remove empty fields
odv2(:,rm_fields) = [];
cols2(rm_fields) = [];

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
  odv2 = movevars(odv2,'date','After','datetime');
  odv2 = movevars(odv2,'time','After','date');
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


%% Parse size bin limits from fieldnames
% split into categories #/L & biovolume
idx_nperL  = find(contains(cols2,'[#/L]') & ~contains(cols2,'='));
idx_biovol = find(contains(cols2,'[mm^3/L]') & ~contains(cols2,'='));
% check that there are equal variables for biovolume, catch and stop if not
if ~isequal(numel(idx_nperL),numel(idx_biovol))
  fprintf('something wrong here..\n')
  keyboard
end
% Pull out size class limits and theoretical bin centers
num_sizes = numel(idx_nperL);
vsize = struct();
vsize.strnam = cell(num_sizes,1);   % string containing only bin limits
vsize.fldnam = cell(num_sizes,1);   % simple field name
vsize.sizemm = nan(num_sizes,2);    % bin range in mm
vsize.sizeum = nan(num_sizes,2);    % bin range in microns
vsize.sizeav_mm = nan(num_sizes,1); % bin center in mm
vsize.sizeav_um = nan(num_sizes,1); % bin center in microns
uvp_diameter_bins_mm = []; % pull out all sizes to check that uvp_diameter_bins correct
for in = 1:num_sizes
  % isolate string containing only bin limits
  vsize.strnam{in} = erase(cols2{idx_nperL(in)},{'(', ')','LPM_','[#/L]','_'});
  % generate basic field to use later
  vsize.fldnam{in} = ['LPM' num2str(in)];
  % pull out sizes from variable name
  if contains(vsize.strnam{in} ,'>') % > 26mm
    str_size1 = erase(vsize.strnam{in},{'mm', '>'});
    sz1 = str2double(str_size1);
    vsize.sizemm(in,:) = [sz1 NaN];
  else
    str_sizes = strsplit(vsize.strnam{in},'-');
    str_size1 = erase(str_sizes{1},{'mm', 'µm'});
    str_size2 = erase(str_sizes{2},{'mm', 'µm'});
    sz1 = str2double(str_size1);
    sz2 = str2double(str_size2);
    % convert numbers to mm
    if contains(vsize.strnam{in},'µm') && contains(vsize.strnam{in},'mm') % only first number in [µm]
      vsize.sizemm(in,:) = [sz1./1000 sz2];
    elseif contains(vsize.strnam{in},'µm') % both numbers in [µm]
      vsize.sizemm(in,:) = [sz1./1000 sz2./1000];
    elseif contains(vsize.strnam{in},'mm') % both numbers in [mm]
      vsize.sizemm(in,:) = [sz1 sz2];
    else
      fprintf('unexpected scenario.. stopping to check here\n')
      keyboard
    end
  end
  % store original sizes
  uvp_diameter_bins_mm = [uvp_diameter_bins_mm vsize.sizemm(in,:)];
  % Calculate size rnage in microns
  vsize.sizeum(in,:) = vsize.sizemm(in,:).*1000; % convert mm to microns
  % Calculate bin center
  % this is really just a theoretical bin center
  % Pull out log10 center of the bin limits in mm
  meansz_mm = logspace(log10(vsize.sizemm(in,1)),log10(vsize.sizemm(in,2)),3); % returns limit1, center, limit2
  vsize.sizeav_mm(in) = round(meansz_mm(2),4); % pull out center of size range
  % Pull out log10 center of the bin limits in microns
  meansz_um = logspace(log10(vsize.sizeum(in,1)),log10(vsize.sizeum(in,2)),3); % returns limit1, center, limit2
  vsize.sizeav_um(in) = round(meansz_um(2),4); % pull out center of size range
  %fprintf(' %s=%s from %.4f-%.4fmm, average size %.4f\n',vsize.fldnam{in},vsize.strnam{in},vsize.sizemm(in,1),vsize.sizemm(in,2),vsize.sizeav(in))
end
% Pull out unique and finize size bin limits
uvp_diameter_bins_mm = unique(uvp_diameter_bins_mm(isfinite(uvp_diameter_bins_mm)));
uvp_diameter_bins_mm = round(uvp_diameter_bins_mm,4);


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
idx_abun = contains(cols2,'[#/L]');
cols2(idx_abun) = strrep(cols2(idx_abun),'LPM_','abundance_'); % strcat('abun_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
cols2(idx_abun) = erase(cols2(idx_abun),'[#/L]');           % erase units from variable name
units(idx_abun) = {'number/L'};

% Pull out biovolume data [#/m^3]
idx_biovol = contains(cols2,'[mm^3/L]');
cols2(idx_biovol) = strrep(cols2(idx_biovol),'LPM_biovolume_','biovolume_'); % strcat('abun_',taxa.ID');  % Build variables as abun_1id,abun_2id... etc
cols2(idx_biovol) = erase(cols2(idx_biovol),'[mm^3/L]');  % erase units from variable name
units(idx_biovol) = {'mm^3/L'};

% add descriptions for variable names 
% First, pull out size bin limits
size_bin_limits = erase(cols2(idx_abun),{'abundance_(' ')' '[#/L]'});  % Pull out text with limits
size_bin_limits = strrep(size_bin_limits,'_µm', ' micrometers'); 
size_bin_limits = strrep(size_bin_limits,'_mm', ' millimeters');
abun_desc   = 'Abundance of particles with an equivalent spherical diameter between';
biovol_desc = 'Biovolume of particles with an equivalent spherical diameter between'; % 'Biovolume, or volume size distribution, for particles';
descr(idx_abun)   = strcat(abun_desc,  {' '},size_bin_limits);
descr(idx_biovol) = strcat(biovol_desc,{' '},size_bin_limits);

% Replace variable descriptions "between >"
descr = strrep(descr,'between >', 'greater than ');

% Replace all µm with um to avoid special character issues
cols2 = strrep(cols2,'µm','um');

% Biovolume fields with no value are exported from ecotaxa as 0, instead of
% blank. That means that when no value flag of -9999 is assigned before
% writing to SeaBASS format, those 0 values stay the same which may lead to
% confusion. Therefore, replace all biovolume == 0 with NaN
biovol = table2array(odv2(:,idx_biovol));
biovol(biovol == 0) = NaN;
odv2(:,idx_biovol)  = array2table(biovol);

% Combine all units into a single string
colunits = strjoin(units,','); %[colstr,cols2{i},','];
% Combine column names into a single string
colstr = strjoin(cols2,','); % colstr=[colstr,cols2{i},','];

% Remove bad data (some profiles seemed to only have a few liters of sampled
% volume per depth bin, resulting in mostly empty particle size
% distributions)
try
  odv2(odv2.LPM_102_128_m___L_1_ == 0,:) = [];
catch
  % may not work in variable names are different
end

%% Generate table that includes all variables with description, units, size bin, etc
fields = table();
fields.VariableName = cols2';
fields.Units = units';
fields.Description = descr';
%% Write ParameterDescriptions table
fprintf('Writing parameter descriptions to %s\n',fullfile(submit_dir,'UVP_PAR_ParameterDescriptions.xlsx'))
writetable(fields,fullfile(submit_dir,'UVP_PAR_ParameterDescriptions.xlsx'))

%% Generate header text to write to file
% first pull out max/min latitude, longitude, date, and time
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
hdr_text = [hdr_text; hdr.comments.PAR; '#end_header'];

%% Specify file writing format
if ~iscell(odv2.site)
  odv2.site = num2str(odv2.site);
end

col_n = length(cols2);
col_n_text = find(contains(cols2,'lat')) - 1; %Latitude is the first DOUBLE column,
fmt_txt = repmat( '%s,', [1 col_n_text] );
fmt_dbl = repmat( '%f,', [1 (col_n - col_n_text)]);
% replace delimiter with new line at end of format string
fmt_dbl(end) = []; fmt_dbl = [fmt_dbl '\n'];
fmt = [fmt_txt fmt_dbl];
% reformat odv2 table to temporary variable to enable simply write to file
odv_write = table2cell(odv2); % convert to cell array to handle different variable types
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
  wfile_name = fullfile(submit_dir,[hdr.ecotaxa_name '_BCO-DMO_PAR.txt']);
catch
  wfile_name = fullfile(submit_dir,strrep(odv_rfile,'.txt','_BCO-DMO_PAR.txt'));
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

end %% MAIN FUNCTION WRITE_BCO_DMO_UVP_PAR








