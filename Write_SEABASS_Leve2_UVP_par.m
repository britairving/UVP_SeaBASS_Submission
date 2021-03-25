function Write_SEABASS_Level2_UVP_par
%% FUNCTION WRITE_SEABASS_LEVEL2_UVP_PAR
% Description:
%   Reads in detailed ODV formatted file exported from Ecotaxa Particle
%   Module and creates a file for submission to NASA's SEABASS database
%   https://seabass.gsfc.nasa.gov/wiki/Data_Submission/
%
% Steps:
%   1. Read in detailed ODV PAR file
%   2. Format fields and column names following SeaBASS format
%   3. Define metadata header that will be printed to .sb file
%   4. Write sb file with DNSD and DVSD
%   5. Write suplementary sb file with non-differential NSD and VSD
%
% References:
%  Picheral M, Colin S, Irisson J-O (2017). EcoTaxa, a tool for the
%  taxonomic classification of images. http://ecotaxa.obs-vlfr.fr.
%
%  UVP_Particle_and_Zooplankton_Characterization_Protocol_20210114.pdf
%
%  Ocean Optics Web Book
%  https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/level-3/creating-particle-size-distributions-data
%
% Authors:
%  Brita Irving <bkirving@alaska.edu>
%  Andrew McDonnell <amcdonnell@alaska.edu>
%% Define script options
include_uncertainty = 0; % 1 = writes uncertainty estimate to file, 0 = does not include uncertainty variables, but does include text on how to include it.
%% Define read and write filenames
%% EXPORTSNP survey cruise on the Sally Ride
cruiseid = 'SR1812'; 
odv_rfile = 'D:\MATLAB\EXPORTS_UVP_SEABASS\export_detailed_20210219_19_10\export_detailed_20210219_19_10_PAR_odv.txt';
odv_wfile = 'EXPORTS-EXPORTSNP_UVP5-ParticulateLevel2_differential_survey_20180814-20180909_R1.sb';
r2r_elog  = 'D:\MATLAB\EXPORTS_UVP_SEABASS\R2R_ELOG_SR1812_FINAL_EVENTLOG_20180913_022931.xlsx';
%% EXPORTSNP process cruise on the Roger Revelle
% cruiseid = 'RR1813'; 
% odv_rfile = 'D:\MATLAB\EXPORTS_UVP_SEABASS\export_detailed_20210128_08_11\export_detailed_20210128_08_11_PAR_odv.txt';
% r2r_elog  = 'D:\MATLAB\EXPORTS_UVP_SEABASS\R2R_ELOG_RR1813_FINAL_EVENTLOG_20180912_170812.xlsx';
% odv_wfile = 'EXPORTS-EXPORTSNP_UVP5-ParticulateLevel2_differential_process_20180814-20180909_R1.sb';

%% Specify write filename
% Release number R0 suggested if data_type=preliminary
% specify write file, this is not computer specific because uses same path
% Write files to submit folder
if ~isdir('submit')
  try
    mkdir('submit');
  catch
    error('cannot create "submit" directory, maybe a permissions error?')
  end
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
% define additional file to store non-differential NSD and VSD
odv_wfile_supp = strrep(odv_wfile,'_differential_','_nondifferential_');

%% Define instrument specs from Calibration documents
% May use these for uncertainty calculation, but not yet set up.
% inst = struct();
% switch cruiseid
%   case 'SR1812'
%     % Values from uvp_calibration_report_of_uvp5sn207a_from_uvp5sn002_20180129.pdf
%     inst.pixel_um  = 0.097;    % Pixel [um]
%     inst.pixel_um2 = 0.009409; % Pixel Area [um^2]
%     inst.im_vol_L  = 1.06;     %Image Volume [L]
%   case 'RR1813'
%     % Values from uvp_calibration_report_of_uvp5sn201a_from_uvp5sn002_20180122_rep_20200325.pdf
%     inst.pixel_um  = 0.094;    % Pixel [um]
%     inst.pixel_um2 = 0.008836; % Pixel Area [um^2]
%     inst.im_vol_L  = 1.13;     %Image Volume [L]
%   otherwise
%     fprintf('pixel area for cruiseid not defined\n')
%     keyboard
% end

%% Read ODV particle data
fprintf('\n  Reading ODV Particle data from... %s\n',odv_rfile)

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
  cols = strrep(cols,')_[mm3_l-1]',')[ppm]');
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

%% Identify all empty variables - i.e. no real values
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

%% DO NOT DELETE 50.8-64_µm and higher size classes
% The 50.8-64_µm size bin will be kept, even if no real values in the field.
% Detection limit for UVP is approximately 100 µm but there is
% commonly data in the 50.8-64µm size bin and we want to maintain all
% columns with data for posterity sake.
% This is also useful because otherwise the size bin limits may be
% misleading, aka metadata header /PSD_bin_size_boundaries.
keep_numL1  = find(contains(cols2,'50.8-64') & contains(cols2,['#/L']));
keep_numL   = [keep_numL1 : find(contains(cols2,['#/L']),1,'last')];
keep_ppm1   = find(contains(cols2,'50.8-64') & contains(cols2,['ppm']));
keep_ppm    = [keep_ppm1 : find(contains(cols2,'[ppm]'),1,'last')];
keep_fields = [keep_numL keep_ppm];
ignore_empty = ismember(rm_fields,keep_fields);
rm_fields(ignore_empty) = [];
%% Now that 50.8-64um are going to be kept, remove empty fields
odv2(:,rm_fields) = [];
cols2(rm_fields) = [];

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
% kept method for R2018a for posterity sake
try % movevars introduced in R2018a
  odv2=movevars(odv2,'date','After','datetime');
  odv2=movevars(odv2,'time','After','date');
catch % if older version of MATLAB
  odv2 = odv2(:,new_order);
end
% sort by datetime then remove the datetime column
odv2.datetime = [];
cols2(i_datetime) = [];
% delete redundant variables that result from reordering
odv2.date_1 = [];
odv2.time_1 = [];
cols2(end-1:end) = []; % original date and time variables at the end


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

%% Parse size bin limits from fieldnames
% split into categories #/L & biovolume
nperL  = find(contains(cols2,'[#/L]') & ~contains(cols2,'='));
biovol = find(contains(cols2,'[ppm]') & ~contains(cols2,'='));
% check that there are equal variables for biovolume, catch and stop if not
if ~isequal(numel(nperL),numel(biovol))
  fprintf('something wrong here..\n')
  keyboard
end
% Pull out size class limits and theoretical bin centers
num_sizes = numel(nperL);
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
  vsize.strnam{in} = erase(cols2{nperL(in)},{'(', ')','LPM_','[#/L]','_'});
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
uvp_diameter_bins_um = uvp_diameter_bins_mm.*1000;

%% Calculate DNSD (Differential number size distribution) [#/m^3/µm] & DVSD (Differential volume size distribution) [µL/m^3/µm]
% NSD = Same as ecotaxa LPM_(size1-size2)[#/L] variables
% VSD = same as ecotaxa LPM_biovolume_(size1-size2)[ppm] variables
% Pull out NSD & VSD & convert to array
NSD = table2array(odv2(:,nperL));  % [#/L]
VSD = table2array(odv2(:,biovol)); % [mm^3/L] or [ppm]
% Biovolume fields with no value are exported from ecotaxa as 0, instead of
% blank. That means that when no value flag of -9999 is assigned before
% writing to SeaBASS format, those 0 values stay the same which may lead to
% confusion. Therefore, replace all biovolume == 0 with NaN
VSD(VSD == 0) = NaN;

% Check if equation for particle volume concentration from Turner et al.,
% 2017 is correct. DOI: 10.1016/j.csr.2017.07.002
% Result: not an exact 1-1 match, I assume this is because the "median
% ESD of each logarithmically-spaced size bin" is calculed differently...
VSD_test = nan(size(VSD));
for nv = 1: size(NSD,2)
  VSD_test(:,nv) = NSD(:,nv) * 4/3 * pi * (vsize.sizeav_mm(nv)/2).^3;  % [#/L]
end
% if ~isequaln(VSD(:,2),VSD_test(:,2))
%   figure,plot(VSD(:,2),'*')
%   hold on
%   plot(VSD_test(:,2),'k.')
%   fprintf('VSD calculation not as expected...\n');
%   keyboard
% end

% Normalize size distribution to the width of the bin
fprintf('Normalizing size distribution to the width of the bin (in mm)\n')
fprintf(' i.e. differential size distribution\n')
for i=1:size(NSD,2)
  wid_mm(i) = diff(vsize.sizemm(i,:));
  DNSD(:,i) = NSD(:,i)./wid_mm(i); % [#/L/mm]
  DVSD(:,i) = VSD(:,i)./wid_mm(i); % [mm^3/L/um] or [ppm/mm]
end

%% Convert units
% Convert DNSD [#/L/mm] to DNSD (Differential number size distribution) [#/m^3/µm]
% units specified by EXPORTS/SEABASS
% convert [#/L/mm] to [#/L/um] where 1mm/1000um =1
% DNSD_1 = DNSD./1000; % [#/L/um]
% convert [#/L/um] to [#/m^3/um] where 1L = 0.001m^3
% DNSD = DNSD_1./0.001; % [#/L/um]
% I.e... this cancels out because  #/L/mm is eqiuvalent to # m^-3 um^-1

% Convert VSDn [mm^3/L/mm] to DVSD (Differential volume size distribution) [µL/m^3/µm]
% 1mm^3 = 1µL
% 1L = 0.001m^3
% 1mm = 1000µm
% I.e... mm^3/L/mm is eqiuvalent to uL m^-3 um^-1
% DVSD = DVSD; % [uL m^-3 um^-1]

%% Calculate Uncertainty
% See https://www.oceanopticsbook.info/view/optical-constituents-of-the-ocean/level-3/creating-particle-size-distributions-data
% For example based, for the large bins, on counting statistics (e.g.
% \sqrt(N)) and at the small end, based also on uncertainties in size. A
% single pixel can include particles smaller than a single pixel but that
% scatter enough to register. Thus, a small particle sitting in between
% pixels can be significantly overestimated in size. This over-estimate
% (+1pixel in every direction) gets relatively smaller the larger the
% particle. Thresholding also induces uncertainties in size, as it defines
% when the signal is large enough to register.
% Since the smallest size bin we are considering is 50.8-64µm and a pixel
% is ~0.01µm, ignore the contribution to uncertainty based on size.

%% Calculate the relative uncertainty 
% relative uncertainty due to counting statistics
rel_unc = sqrt( NSD .* odv2.SampledVolume_L_) ./ (NSD .* odv2.SampledVolume_L_); % [%]
u_NSD  = rel_unc .* NSD;  % [%]*[#/L]       = [#/L]    convert to [#/m^3] later
u_VSD  = rel_unc .* VSD;  % [%]*[mm^3/L]    = [mm^3/L] convert to [µL/m^3] later
u_DNSD = rel_unc .* DNSD; % [%]*[#/m^3/µm]  = [#/m^3/µm]
u_DVSD = rel_unc .* DVSD; % [%]*[µL/m^3/µm] = [µL/m^3/µm]

rel_unc_text = {'!  relative uncertainty [none]';...
                '!  relative_unc = sqrt(PSD_NSD_{#}umsize*volume)/(PSD_NSD_{#}umsize*volume)'};
u_NSD_text   = {'!  uncertainty in the number size distribution [#/m^3]';...
                '!  PSD_NSD_{#}umsize_unc = relative_unc*PSD_NSD_{#}umsize'};
u_DNSD_text  = {'!  uncertainty in the differential number size distribution [#/m^3/µm]';...
                '!  PSD_DNSD_{#}umsize_unc = relative_unc*PSD_DNSD_{#}umsize'};
u_VSD_text   = {'!  uncertainty in the volume size distribution [µL/m^3]';...
                '!  PSD_VSD_{#}umsize_unc = relative_unc*PSD_VSD_{#}umsize'};
u_DVSD_text  = {'!  uncertainty in the differential volume size distribution [µL/m^3/µm]';...
                '!  PSD_DVSD_{#}umsize_unc = relative_unc*PSD_DVSD_{#}umsize'};
u_text_differential    = [rel_unc_text; u_DNSD_text; u_DVSD_text];
u_text_nondifferential = [rel_unc_text; u_NSD_text;  u_VSD_text];
%% Calculate uncertainty based on parameter derivations
% % Based on the above advice from Emmanuel Boss, first calculate uncertainty
% % simply based on counting statistics
% % Multply number concentration by sample volume to get #
% N   = NSD .* odv2.SampledVolume_L_;
% u_N = sqrt( N );
%
% % Now calculate differential fields using uncertainty to get same units
% u_NSD = u_N ./ odv2.SampledVolume_L_;  % [#/L]
% 
% % Normalize uncertainty of size distribution to the width of the bin
% for ns=1:size(u_NSD,2)
%   u_DNSD(:,ns) = u_NSD(:,ns)./wid_mm(ns); % [#/L/mm]
% end
% 
% fprintf('How to calculate uncertainty in volume estimatation... getting back to uncertainty in size...\n')
% %Turner, J.S., J.L. Pretty, and A.M.P. McDonnell. (2017). "Marine particles
% %in the Gulf of Alaska shelf system: spatial patterns and size
% %distributions from in situ optics". Continental Shelf Research. 145:13–20.
% % doi: https://doi.org/10.1016/j.csr.2017.07.002
% % " The UVP generated particle size distributions in terms of a count per
% %   volume of water photographed (#/l). These values were then converted to a
% %   volumetric concentration (µl/l) using the median ESD (mm) of each
% %   logarithmically-spaced size bin and assuming spherical shape using the
% %   following equation"
% %   V = N*4/3*pi*(d/2)^3
% %   For each UVP size class, Vis the particle volume concentration (µl/ l),
% %   N is the number of particles counted in each liter of water (#/l), and
% %   d is the median diameter of the respective size bin.
% % Now calculate differential fields using uncertainty to get same units%
% u_VSD = nan(size(VSD));
% for ns = 1:size(NSD,2)
%   % calculate uncertainty in volume particle concentration
%   u_VSD(:,ns) = u_NSD(:,ns) * 4/3 * pi * (vsize.sizeav_mm(ns)/2).^3;  % [mm^3/L]
%   % uncertainty in differential field
%   u_DVSD(:,ns) = u_VSD(:,ns) ./ wid_mm(ns); % [mm^3/L/mm] == [µL/m^3/µm]
% end
% u_NSD_text  = {'!  uncertainty in the number size distribution';...
%                '!  PSD_NSD_{PSD_bin_size}_unc = sqrt(PSD_NSD_{PSD_bin_size}*volume)/volume'};
% u_DNSD_text = {'!  uncertainty in the differential number size distribution';...
%                '!  PSD_DNSD_{PSD_bin_size}_unc = ( sqrt(PSD_NSD_{PSD_bin_size}*volume)/volume ) / ( PSD_bin_size_boundaries2-PSD_bin_size_boundaries1 )'};
% u_VSD_text  = {'!  uncertainty in the volume size distribution';...
%                '!  PSD_VSD_{PSD_bin_size}_unc  = ( sqrt(PSD_NSD_{PSD_bin_size}*volume)/volume ) * 4/3 * pi * (PSD_bin_size/2)^3'};
% u_DVSD_text = {'!  uncertainty in the differential volume size distribution';...
%                '!  PSD_DVSD_{PSD_bin_size}_unc = ( ( sqrt(PSD_NSD_{PSD_bin_size}*volume)/volume ) * 4/3 * pi * (PSD_bin_size/2)^3 ) / ( PSD_bin_size_boundaries2-PSD_bin_size_boundaries1 )'};
% u_text_differential    = ['!  PSD_bin_size = Median bin diameter in microns'; u_DNSD_text; u_DVSD_text];
% u_text_nondifferential = ['!  PSD_bin_size = Median bin diameter in microns'; u_NSD_text;  u_VSD_text];

%% Add PSD_DNSD [number/m^3/um] to data table & update column name
% Differential number size distribution. Must be used in conjunction with a
% suffix to specify sizes, e.g., PSD_DNSD_###umsize where ### in
% "_###umsize" = Median bin diameter in microns.
% Units are added later
num_size_classes = size(DNSD,2)-1; % Do not include last size class because >26mm so differential is NaN
for i=1:num_size_classes
  % Generate new variable name following SeaBASS format
  % since size is specified by median bin diameter, which is really a
  % theoretical bin center, just round the log10 center of the size class
  nsd_name = ['PSD_DNSD_' num2str(round(vsize.sizeav_um(i))) 'umsize'];
  % add to existing column name cell
  cols2 = [cols2, nsd_name];
  % add variable to data table
  odv2.(['DNSD' num2str(i)]) = DNSD(:,i);
  if include_uncertainty
    cols2 = [cols2, [nsd_name '_unc']];
    % add uncertainty variable to data table
    odv2.(['u_DNSD' num2str(i)]) = u_DNSD(:,i);
  end
end

%% Add PSD_DVSD [uL/m^3/um] to data table & update column name
% Differential volume size distribution. Must be used in conjunction with a
% suffix to specify sizes, e.g., PSD_DVSD_###umsize where ### is the Median
% bin diameter in microns.
% Units are added later
for i=1:num_size_classes
  % Generate new variable name following SeaBASS format
  % since size is specified by median bin diameter, which is really a
  % theoretical bin center, just round the log10 center of the size class
  vsd_name = ['PSD_DVSD_' num2str(round(vsize.sizeav_um(i))) 'umsize'];
  % add to existing column name cell
  cols2 = [cols2, vsd_name];
  % add variable to data table
  odv2.(['DVSD' num2str(i)]) = DVSD(:,i);
  if include_uncertainty
    cols2 = [cols2, [vsd_name '_unc']];
    % add uncertainty to data table
    odv2.(['u_DVSD' num2str(i)]) = u_DVSD(:,i);
  end
end

%% Remove original LPM variables
% initialize new variables so don't overwrite previous
odv3  = odv2;
cols3 = cols2;
% find variables that we want to remove
remthis = find(contains(cols3,'LPM'));
odv3(:,remthis) = [];
cols3(remthis)  = [];

%% Format SeaBASS variable name and units strings
% format based on SeaBASS requirements
% https://seabass.gsfc.nasa.gov/wiki/stdfields
units = cell(size(cols3));
for i=1:length(cols3)
  % BUILD UNITS
  k1 = strfind(cols3{i},'[');
  k2 = strfind(cols3{i},']');
  if contains(cols3{i},'PSD_DVSD')    % differential volume size distribution
    % this will also account for the uncertainties
    units{i} = 'uL/m^3/um';
  elseif contains(cols3{i},'PSD_DNSD')% differential number size distribution
    % this will also account for the uncertainties
    units{i} = 'number/m^3/um';
  elseif contains(cols3{i},{'lat' 'lon'},'IgnoreCase',true) % latitude and longitude
    units{i} = 'degrees';
    % remove units string from column name
    cols3{i}=cols3{i}(1:k1-2);
  elseif contains(cols3{i},'Depth','IgnoreCase',true) % depth
    units{i} = 'm';
    % remove units string from column name
    cols3{i}=cols3{i}(1:k1-2);
  elseif isempty(k1) % Pull out other cases
    if strcmp(cols3{i},'date')
      units{i} = 'yyyymmdd';
    elseif strcmp(cols3{i},'time')
      units{i} = 'hh:mm:ss';
    elseif strcmp(cols3{i},'R2R_Event')
      units{i} = 'none';
    else
      units{i} ='none';
    end
  else
    fprintf('  %s check units are correct: %s\n',cols3{i},cols3{i}(k1+1:k2-1))
    %fprintf('In pause mode... press any key to continue\n')
    %pause
    %extract only the units string
    units{i} = cols3{i}(k1+1:k2-1);
    % remove units string from column name
    cols3{i}=cols3{i}(1:k1-2);
  end
  
  % FORMAT COLUMN NAMES
  % remove extraneous text from strings containing variable names
  cols3{i}=strrep(cols3{i},'_[*]','');
  %Correct some variable names to match with SEABASS requirements
  cols3{i}=strrep(cols3{i},'Site','station');
  cols3{i}=strrep(cols3{i},'Station','station_alt_id');
  cols3{i}=strrep(cols3{i},'Latitude','lat');
  cols3{i}=strrep(cols3{i},'Longitude','lon');
  % Variable SHOULD be bin_depth but fcheck throws and error
  % "Measurement depth required as either a header value or a column in the
  %  data block.  All data MUST have an associated depth"
  %cols3{i}= strrep(cols3{i},'Depth','bin_depth'); % bin_depth = Nominal or center depth for each data bin
  cols3{i}= strrep(cols3{i},'Depth','depth');  % depth = Depth of measurement
  cols3{i}=strrep(cols3{i},'Sampled_volume','volume');
end
% Combine all units into a single string
colunits = strjoin(units,','); %[colstr,cols3{i},','];
% Combine column names into a single string
colstr = strjoin(cols3,','); % colstr=[colstr,cols3{i},','];

% Correct some unit names to match with SEABASS requirements
colunits=strrep(colunits,'degrees_north','degrees');
colunits=strrep(colunits,'degrees_east','degrees');

% Remove bad data (some profiles seemed to only have a few liters of sampled
% volume per depth bin, resulting in mostly empty particle size
% distributions)
odv3(odv3.DNSD5==0,:)=[];

%% Generate SeaBASS header text
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

% First - grab structure with basic metadata from SeaBASS_metadata
[~,original_file,ext] = fileparts(odv_rfile);
original_file = [original_file ext];
hdr = SeaBASS_metadata_headers(cruiseid,original_file);

% Write formatted header for sb file
hdr_SEABASS = format_UVP_SeaBASS_metadata_header;

%% Specify file writing format
col_n = length(cols3);
col_n_text = find(contains(cols3,'lat')) - 1; %Latitude is the first DOUBLE column,
fmt_txt = repmat( '%s,', [1 col_n_text] );
fmt_dbl = repmat( '%f,', [1 (col_n - col_n_text)]);
% replace delimiter with new line at end of format string
fmt_dbl(end) = []; fmt_dbl = [fmt_dbl '\n'];
fmt = [fmt_txt fmt_dbl];
% reformat odv2 table to temporary variable to enable simply write to file
odv_write = table2cell(odv3); % convert to cell array to handle different variable types
odv_write = odv_write';            % transpose because fprintf function prints data columnwise
% convert NaNs to missing identifier
% odv_write(cellfun(@(x) isnumeric(x) && isnan(x), odv_write)) = {''}; % convert NaNs to '' for ODV format
odv_write(cellfun(@(x) isnumeric(x) && isnan(x), odv_write)) = {-9999}; % missing=-9999 SeaBASS required header field

%% Write odv table to DIFFERENTIAL file
fprintf('\n  Deleting data file, if it exist (if it does not exist, you may get a WARNING, but this is OK: %s\n',odv_wfile);
eval(['delete ' odv_wfile])

fprintf('\n  Writing data in ODV format to: %s\n',fullfile('submit',odv_wfile));
fileID = fopen(fullfile('submit',odv_wfile),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile('submit',odv_wfile))
  keyboard
end
fprintf(fileID,'%s\n',hdr_SEABASS{:});  % write header
fprintf(fileID,fmt,odv_write{:});     % write data
fclose(fileID);                 % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
%fprintf('%s\n',hdr_SEABASS{:}) % write header
%fprintf(fmt,odv_write{:})      % write data

%% Write file with non-differential NSD and VSD
% Write additional file that includes NSD and VSD. 
% Units of differential fields
% PSD_DNSD [#/m^3/µm]
% PSD_DVSD [µL/m^3/µm]
% Therefore.. convert 
% NSD from [#/L]    to [#/m^3]  --> NSD * 1000.0;
% VSD from [mm^3/L] to [µL/m^3] --> VSD * 1000.0; (1uL=1mm^3, 1000L=1m^3)

odv4 = odv3;
cols4 = cols3;

% get rid of DNSD and DVSD
rm_after = find(contains(cols4,'PSD_'),1);
odv4(:,rm_after:end) = [];
cols4(rm_after:end)  = [];
units(rm_after:end)  = [];

% Add NSD and VSD
% SeaBASS has PSD as field, so do not use PSD and instead use NSD
% ---------------------------------------------------------------
% Add NSD [number/L] to data table & update column name
for i=1:num_size_classes
  % Generate new variable name following SeaBASS format
  % since size is specified by median bin diameter, which is really a
  % theoretical bin center, just round the log10 center of the size class
  nsd_name = ['PSD_NSD_'  num2str(round(vsize.sizeav_um(i))) 'umsize'];
  % add to existing column name cell
  cols4 = [cols4 nsd_name];
  % add unit to units structure
  %units = [units 'number/L'];
  % convert from number/L to number/m^3
  units = [units 'number/m^3'];
  % add variable to data table
  odv4.(['NSD' num2str(i)]) = NSD(:,i).*1000.0; % 1000L/m^3
  if include_uncertainty
    % add name and unit
    cols4 = [cols4 [nsd_name '_unc']];
    units = [units  'number/m^3'];
    % add variable to data table
    % u_VSD was in units of [mm^3/L], so convert to [uL/m^3]
    odv4.(['u_NSD' num2str(i)]) = u_NSD(:,i).*1000.0; % 1uL=1mm^3, 1000L=1m^3
  end
end
% ---------------------------------------------------------------
% Add VSD [mm^3/L] to data table & update column name
% convert mm^3/L to uL/m^3
for i=1:num_size_classes
  % Generate new variable name following SeaBASS format
  % since size is specified by median bin diameter, which is really a
  % theoretical bin center, just round the log10 center of the size class
  vsd_name = ['PSD_VSD_'  num2str(round(vsize.sizeav_um(i))) 'umsize'];
  % add to existing column name cell
  cols4 = [cols4 vsd_name];
  % add unit to units structure
  %units = [units 'mm^3/L'];
  units = [units 'uL/m^3'];
  % add variable to data table
  odv4.(['VSD' num2str(i)]) = VSD(:,i).*1000.0; % 1uL=1mm^3, 1000L=1m^3
  if include_uncertainty
    % add name and unit
    cols4 = [cols4 [vsd_name '_unc']];
    units = [units 'uL/m^3'];
    % add variable to data table
    % u_VSD was in units of [mm^3/L], so convert to [uL/m^3]
    odv4.(['u_VSD' num2str(i)]) = u_VSD(:,i).*1000.0; % 1uL=1mm^3, 1000L=1m^3
  end
end

% Combine all units into a single string
colunits = strjoin(units,','); %[colstr,cols3{i},','];
% Combine column names into a single string
colstr = strjoin(cols4,','); % colstr=[colstr,cols3{i},','];

%% Replace units and fields in header
hdr_SEABASS(contains(hdr_SEABASS,'/fields')) = {['/fields=' colstr]};
hdr_SEABASS(contains(hdr_SEABASS,'/units')) = {['/units=' colunits]};

%% Swap 'data_file_name' and 'associated_files' in metadata header
hdr_SEABASS(contains(hdr_SEABASS,['/data_file_name=' odv_wfile])) = {['/data_file_name=' odv_wfile_supp]};
hdr_SEABASS(contains(hdr_SEABASS,['/associated_files=' odv_wfile_supp])) = {['/associated_files=' odv_wfile]};
% Swap uncertainty descriptions
hdr_SEABASS(contains(hdr_SEABASS,u_text_differential)) = u_text_nondifferential;

%% Specify file writing format
col_n = length(cols4);
col_n_text = find(contains(cols4,'lat')) - 1; %Latitude is the first DOUBLE column,
fmt_txt = repmat( '%s,', [1 col_n_text] );
fmt_dbl = repmat( '%f,', [1 (col_n - col_n_text)]);
% replace delimiter with new line at end of format string
fmt_dbl(end) = []; fmt_dbl = [fmt_dbl '\n'];
fmt = [fmt_txt fmt_dbl];
% reformat data table to temporary variable to enable simply write to file
odv_write = table2cell(odv4); % convert to cell array to handle different variable types
odv_write = odv_write';       % transpose because fprintf function prints data columnwise
% convert NaNs to missing identifier
odv_write(cellfun(@(x) isnumeric(x) && isnan(x), odv_write)) = {-9999}; % missing=-9999 SeaBASS required header field

%% Write odv table to NONDIFFERENTIAL file
fprintf('\n  Writing data in ODV format to: %s\n',fullfile('submit',odv_wfile_supp));
fileID = fopen(fullfile('submit',odv_wfile_supp),'w');  % open file
if fileID < 0
  fprintf(' *** Error opening file %s\n',fullfile('submit',odv_wfile_supp))
  keyboard
end
fprintf(fileID,'%s\n',hdr_SEABASS{:}); % write header
fprintf(fileID,fmt,odv_write{:});      % write data
fclose(fileID);                 % close file
% To troubleshoot why not printing correctly, comment above and just print to screen
%fprintf('%s\n',hdr_SEABASS{:}) % write header
%fprintf(fmt,odv_write{:})      % write data

fprintf('workflow is finished\n')

%% FUNCTION META = FORMAT_UVP_SEABASS_METADATA_HEADER
  function hdr_SEABASS = format_UVP_SeaBASS_metadata_header
    hdr_SEABASS={'/begin_header';
      ['/investigators='  hdr.investigators];
      ['/affiliations='   hdr.affiliations];
      ['/contact='        hdr.contact];
      ['/experiment='     hdr.experiment];
      ['/cruise='         hdr.cruise];
      ['/station='        hdr.station];
      ['/data_file_name=' odv_wfile];
      ['/associated_files=' odv_wfile_supp]; % associated file that contains non-differential NSD and VSD
      ['/original_file_name=' original_file]; % The original name of the data file, if different from the current /data_file_name. Designed to be a reference for the contributor.
      ['/documents='    strjoin(hdr.documents,',')];
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
      ['/missing='                  hdr.missing];
      ['/delimiter='                hdr.delimiter];
      ['/instrument_manufacturer='  hdr.inst_mfr];
      ['/instrument_model='         hdr.inst_model];
      ['/calibration_files='        hdr.calfiles];
      ['/calibration_date='         hdr.caldates];
      '/PSD_bin_size_method=10_to_the_power_of_the_mean_of_the_log10_transformed_size_bin_limits_in_micrometers'; %specify method used to select the nominal bin size such as arithmetic_mean or geometric_mean, just make sure there are no spaces in the string
      '! The sizes reported in the field names represent the bin center, which is equal to 10 to the power of the mean of the individual log10 transformed size bin limits.'
      ['/PSD_bin_size_boundaries=' char(strjoin(string(uvp_diameter_bins_um),','))];  %please provide a comma-separated list with the bin size boundaries in increasing order, e.g., 5,9.5,15,20 ...
      '! size of the bin limits in micrometers';
      '!'};
    if include_uncertainty
        hdr.comments = [hdr.comments; ...
          '!';...
          '! Uncertainties based on counting statistics';...
          u_text_differential];
    else % Include description of how uncertainties could be calculated
        hdr.comments = [hdr.comments; ...
          '!';...
          '! Uncertainties can be calculated based on counting statistics';...
          '! For example:' ;...
          u_text_differential];
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








