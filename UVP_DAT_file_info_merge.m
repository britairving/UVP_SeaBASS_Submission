function UVP_DAT_file_info_merge
%% FUNCTION UVP_DAT_FILE_INFO_MERGE
% Description:
%   Reads DAT files from the UVP "work" folder
%
% Author:
%  Brita Irving <bkirving@alaska.edu>
%%
write_filename = '/Users/bkirving/Documents/MATLAB/UVP_project_data/SR1812/uvp5_sn207_2018_exports_np_sr1812_DAT.csv';

work_directory = '/Volumes/MPDL/UVP_uvpb/uvp5_sn207_2018_exports_np_sr1812/work';
cast_folders = dir(work_directory);
% remove erroneous folders
bad = contains({cast_folders.name},'.');
cast_folders(bad) = [];

%% Define column titles
% from UVP5_User_Manual_20210201 section 2.5.4 "DAT files description"
cols = {'index' 'image' ...
  'p_centibar' 'angle_180' 'angle_1802' 'mb_tempC' 'batt_mV' 'aaaaa2' 'aaaaa3' 'aaaaa4' 'aaaaa5' 'cd_temp' 'camera_temp' 'aaaaa8'... % Sensor data
  'nb_blobs_PG' 'mean_area_PG' 'mean_grey_PG' 'nb_blobs_G' 'mean_grey_G'};

%% Initialize dat structure
dat = struct();
dat.folder = {cast_folders.name};

%% Loop through folders and read DAT info
for ncast = 1:numel(cast_folders)
  castdir = fullfile(work_directory,cast_folders(ncast).name);
  datfile = dir(fullfile(castdir,'*.dat'));
  if numel(datfile) == 0
    fprintf('NO dat files in this folder! %s\n',cast_folders(ncast).name)
  else
    for ndat = 1:numel(datfile)
      try
        filename = fullfile(datfile(ndat).folder,datfile(ndat).name);
        fprintf('Reading file %s\n',filename)
        t = readtable(filename,'FileType','text','Delimiter',';');
        %if size(t,2) > 19
        %  t = t(:,1:19);
        %end
        % add column names
        %t.Properties.VariableNames = cols;
        % ONLY keep the first three columns
        % index, image, and pressure in centibars
        t = t(:,1:3); 
        % add column names
        t.Properties.VariableNames = cols(1:3);
        % add 'datfilename' and 'work_folder' 
        t.datfilename = repmat({datfile(ndat).name},size(t,1),1);
        t.work_folder = repmat({cast_folders(ncast).name},size(t,1),1);
        if ~isfield(dat,'data')
          dat.data = t;
          dat.filename = {datfile(ndat).name};
        else
          dat.data = [dat.data; t];
          dat.filename = [dat.filename; {datfile(ndat).name}];
        end
      catch
        fprintf('could not add \n')
        keyboard
      end % try catch 
    end % Loop through all dat files and read in data
  end % If folder contains a .dat file, read it
end % Loop through all work folders 

fprintf('DONE!\n')
keyboard

%% Write abbreviated DAT information to file 
fprintf('Writing data to %s\n',write_filename)
writetable(dat.data,write_filename,'FileType','text','Delimiter',',');

end %% MAIN FUNCTION