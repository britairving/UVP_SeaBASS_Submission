function odv = UVP_merge_R2R(odv,r2r_elog,cruiseid)
%% FUNCTION UVP_MERGE_R2R
% Description:
%   Finds matching R2R_Event based on cruiseID
%
% Authors:
%  Brita Irving <bkirving@alaska.edu>

%% Read R2R log
  if contains(r2r_elog,'.xlsx')
    r2r = readtable(r2r_elog,'FileType','spreadsheet');
  else
    r2r = readtable(r2r_elog,'FileType','text');
  end
  % Initialize r2r event field
  odv.R2R_Event = cell(size(odv.cruise));
  odv.R2R_Station = cell(size(odv.cruise));
  %odv.r2r_cast  = cell(size(odv.cruise));
  
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
    %% SG2105 - EXPORTS North Atlantic SdG Sarmiento de Gamboa cruise
    case 'SG2105'
      % NOTE 'SG2105_UVP6' -- UVP6 LP deployed on float TZEX, R2R_Event     '20210507.0713.001'    2.0211e+07     '07-May-2021 07:12:57'    'TZEX'                           'other'              NaN       'NaN'       'NaN'      48.945      -14.825     ''            'eCeballos-Romero1'    'Surface deploy and recovery test'                                                           

      try
        % store lat and long to check
        R2R_Latitude  = nan(size(odv.cruise));
        R2R_Longitude = nan(size(odv.cruise));
        [sites,iu] = unique(odv.site,'stable');
        
        % remove empty entry (necessary because odv format only has single
        % entry per profile
        if any(strcmp(sites,'') == 1)
          rm_empty = strcmp(sites,'');
          sites(rm_empty) = [];
          iu(rm_empty) = [];
        end
        % pull out profiles from uvp data
        profiles1 = erase(odv.profile(iu),{'sg2105_00' 'sg2105_0' 'sg2105_'});
        % loop through profiles and pull out r2revent
        for nr = 1:numel(sites)
          % math by UVP rawfilename, then station or cast
          mtch = strcmp(r2r.Cast,profiles1{nr}) & strcmp(r2r.Instrument,'CTD911');
          if sum(mtch) == 1
            odv.R2R_Event{iu(nr)} = r2r.Event{mtch};
            odv.R2R_Station{iu(nr)} = r2r.Station{mtch};
            % replace UVP latitude/longitude with R2R lat lon
            odv.latitude(iu(nr))  = r2r.Latitude(mtch);
            odv.longitude(iu(nr)) = r2r.Longitude(mtch);
          elseif sum(mtch) >= 2
            mtch = strcmp(r2r.Cast,profiles1{nr}) & strcmp(r2r.Instrument,'CTD911') & strcmp(r2r.Action,'deploy');
            if sum(mtch) > 1
              mtch = find(mtch == 1,1,'last');
            end
            odv.R2R_Event{iu(nr)}   = r2r.Event{mtch};
            odv.R2R_Station{iu(nr)} = r2r.Station{mtch};
            % replace UVP latitude/longitude with R2R lat lon
            odv.latitude(iu(nr))  = r2r.Latitude(mtch);
            odv.longitude(iu(nr)) = r2r.Longitude(mtch);
          else
            fprintf('r2r event not found, stopping here\n')
            keyboard
          end
        end
      catch
        keyboard
      end
      % Remove R2R_Station 
      odv.R2R_Station = [];
    %% DY131 - EXPORTS North Atlantic Survey ship RRS Discovery  
    case 'DY131'
      try
      % store lat and long to check
      R2R_Latitude  = nan(size(odv.cruise));
      R2R_Longitude = nan(size(odv.cruise));
      [sites,iu] = unique(odv.site,'stable');
      
      % remove empty entry (necessary because odv format only has single
      % entry per profile
      if any(strcmp(sites,'') == 1)
        rm_empty = strcmp(sites,'');
        sites(rm_empty) = [];
        iu(rm_empty) = [];
      end
      % pull out profiles from uvp data
      profiles1 = erase(odv.profile(iu),{'ctd00' 'ctd0' 'ctd'});
      % loop through profiles and pull out r2revent
      for nr = 1:numel(sites)
        % math by UVP rawfilename, then station or cast
        mtch = strcmp(r2r.R2R_Cast,profiles1{nr}) & strcmp(r2r.R2R_Instrument,'CTD911');
        if sum(mtch) == 1
          odv.R2R_Event{iu(nr)} = r2r.R2R_Event{mtch};
          odv.R2R_Station{iu(nr)} = r2r.R2R_Station{mtch};
          % replace UVP latitude/longitude with R2R lat lon
          odv.latitude(iu(nr))  = r2r.R2R_Latitude(mtch);
          odv.longitude(iu(nr)) = r2r.R2R_Longitude(mtch);
          %R2R_Latitude(iu(nr))  = r2r.R2R_Latitude(mtch);
          %R2R_Longitude(iu(nr)) = r2r.R2R_Longitude(mtch);
        elseif sum(mtch) >= 2 
          mtch = strcmp(r2r.R2R_Cast,profiles1{nr}) & strcmp(r2r.R2R_Instrument,'CTD911') & strcmp(r2r.R2R_Action,'deploy');
          if sum(mtch) > 1
            mtch = find(mtch == 1,1,'last');
          end
          odv.R2R_Event{iu(nr)}   = r2r.R2R_Event{mtch};
          odv.R2R_Station{iu(nr)} = r2r.R2R_Station{mtch};
          % replace UVP latitude/longitude with R2R lat lon
          odv.latitude(iu(nr))  = r2r.R2R_Latitude(mtch);
          odv.longitude(iu(nr)) = r2r.R2R_Longitude(mtch);
          %           R2R_Latitude(iu(nr))  = r2r.R2R_Latitude(mtch);
          %           R2R_Longitude(iu(nr)) = r2r.R2R_Longitude(mtch);
          %           if abs(odv.latitude(iu(nr)) - R2R_Latitude(iu(nr))) > 0.1
          %             fprintf('R2R |%s \t %f, %f\n',r2r.R2R_Date{mtch}, round(R2R_Latitude(iu(nr)),3),round(R2R_Longitude(iu(nr)),3))
          %             fprintf('UVP |%s \t\t %f, %f\n',datestr(odv.datetime(iu(nr)),'ddd  dd mmm yyyy HH:MM:SS'),round(odv.latitude(iu(nr)),3),round(odv.longitude(iu(nr)),3))
          %           end
        else
          fprintf('r2r event not found, stopping here\n')
          keyboard
        end
      end
      catch
        keyboard
      end
    otherwise
      fprintf('R2R event log not match up yet, ignoring for now\n')
      odv.R2R_Event = [];
  end
end