function taxa = WoRMS_AphiaID_taxa_match(taxa,check_children_manual,data_provider)
%% FUNCTION TAXA = WORMS_APHIAID_TAXA_MATCH(TAXA)
% Description:
%   Calls the WoRMS webservice from Matlab to match taxonomic name with
%   WoRMS AphiaID.  
%
% Steps:
%   1. Check for direct match with getAphiaID 
%      getAphiaID Returns:
%         NULL when no match is found
%         -999 when multiple matches are found
%         an integer (AphiaID) when one exact match was found
%   2. If multiple matches are found (-999) & parent name available, it
%      matches the parent, then grandparent. Grandparent is found by
%      searching Names for a match with Name_parent, etc. 
%   3. IF check_children_manual (Default=true)
%      if no match is found (NULL) you can list all children of the parent
%      and choose to select a child 
%      For example:
%         Teleostei_X --> Teleostei 
%         Trachylina  --> Trachylinae
%      Once you select a child - be sure to add it to section "STEP3 |
%      Check known special cases" so it is automatic in the future. 
%
% Inputs:
%   taxa | table with fields Name, Name_parent, & Name_original
%     For Example:
%     taxa = table();
%     taxa.Name          = {}; % taxa name, lowest level
%     taxa.Name_parent   = {}; % taxa parent name, if available
%     taxa.Name_original = {}; % name of original variable
%
%   check_children_manual | true/false boolean
%     true  = checks full classifcation of parent if child doesn't have a
%             match and asks user.
%     false = skips this.
%
%   data_provider | Name of database where taxonomic names were generated.
%     By passing this variable, you can account for common names that are
%     known to not match and AphiaID in WoRMS. 
%     For Example: data_provider = 'ecotaxa'; 
%         Taxonomic names from manual or autmoatic identification in
%         Ecotaxa <https://ecotaxa.obs-vlfr.fr/>. Ecotaxa IDs can often
%         contain parts, such as leg, colonial, spines, etc...
% 
% Outputs:
%   taxa | table with added fields AphiaID, AphiaID_parent, scientificNameID, and Comment.
%   Example
%     taxa.AphiaID          = {}; % AphiaID (NULL = no match)
%     taxa.AphiaID_parent   = {}; % AphiaID of parent (NULL = no match)
%     taxa.scientificNameID = {}; % [Network identifer]:[Reference database]:[Namespace]:[Taxon identifer]
%     taxa.Comment          = {}; % Comments
%
% References:
%  OCB_PTWG_TM_FinalDraft_Dec2020.docx
%
% Author:
%  Brita Irving <bkirving@alaska.edu>
%% Set debug mode
dbstop if error
%% STEP0 | Define options used in script
% check_children_manual
if nargin == 1
  check_children_manual = true; % Default=true
end

% Define arguements needed for WoRMS webservice functions
args = struct();
args.marine_only = true; % marine_only: Limit to marine taxa. Default=true
args.like        = true; % like: Add a '%'-sign after the ScientificName (SQL LIKE function). Default=false
args.fuzzy       = true; % fuzzy: This parameter is deprecated (and ignored since '2013-07-17'). Please use the function matchAphiaRecordsByNames() for fuzzy/near matching
args.offset      = true; % offset: Starting recordnumber, when retrieving next chunk of (50) records. Default=1

% Boolean if parent names are available
if ismember('Name_parent',taxa.Properties.VariableNames) 
  flag_parent_available = true;
else
  flag_parent_available = false;
end

% Define new fields in table
taxa.AphiaID = cell(size(taxa.Name));
taxa.scientificName   = cell(size(taxa.Name));
if flag_parent_available
  taxa.AphiaID_parent = cell(size(taxa.Name));
end
taxa.scientificNameID = cell(size(taxa.Name));
taxa.Comment          = cell(size(taxa.Name));

%% STEP1 | Define common/known issues related to specific data_provider
if nargin == 3
  if ischar(data_provider)
    switch data_provider
      
      case {'ecotaxa' 'Ecotaxa' 'ECOTAXA'} % match all cases if spelling is different
        % Add special known cases
        special = struct();
        special.Teleostei     = 'Teleostei X';    % Likely named as temporary field during manual Ecotaxa identification (EXPORTS SR1812 cruise)
        special.Trachylinae   = 'Trachylina';  % Might be common misspelling in Ecotaxa
        special.Globigerinina = 'Globigerinacea';
        special.Rhizaria      = {'Rhizaria X sp.' 'rhizaria like'};
        special.Ctenophora    = {'Ctenophora_X' 'Ctenophora_XX' 'Ctenophora_sp.'};
        
        % IF the data provider is confident that the ROI is a living particle
        % but cannot be identified to a specific taxonomic rank, then it should
        % be classified to the rank of Eukaryota (see Subchapter 2.4.1 for an
        % example).
        
      otherwise
        fprintf('data_provider has not been set to match parts to the parent yet\n')
    end
  else
    error('Expected data_provider to be a character string, e.g. "ecotaxa"')
  end
end

%% STEP2 | Set up WoRMS webservice and create object
% Follows instructions for WoRMS webservice for Matlab
%   <http://www.marinespecies.org/aphia.php?p=webservice&type=matlab>
%
% Step 1: Download the WSDL file
%    "http://www.marinespecies.org/aphia.php?p=soap&wsdl=1"
%    Name the file and put in folder of choice (e.g. "rootdir/WoRMS/aphia.xml")
aphia_filename = 'WoRMS_Aphia.xml';
if ~exist(aphia_filename,'file')
  try % Instead of copying website directly, try to automate with Matlab
    aphia_filename = websave('WoRMS_Aphia.xml','https://www.marinespecies.org/aphia.php?p=soap&wsdl=1');
  catch
    fprintf('Could not automatically save WSDL file with websave.m builtin function\n')
    fprintf('You will need to do this manually\n')
    fprintf('\n');
    fprintf('Step 1: Download the WSDL file\n')
    fprintf('  Save "http://www.marinespecies.org/aphia.php?p=soap&wsdl=1" to "%s"\n',aphia_filename)
    fprintf('  Entering keyboard mode so you can do this\n');
    keyboard
  end
end
% Step 2: in Matlab: move to folder where wsdl file is located and run:
try
  createClassFromWsdl(aphia_filename);
catch
  fprintf('Could not create Class from Wsdl, check if file exists and is in the path: %s\n',aphia_filename);
  error('Can not continue until class is created so can call object');
end
% Step 4: in Matlab: create an object of the class:
obj = AphiaNameService;

%% STEP3 | Query AphiaID and parent AphiaID
num_taxa = size(taxa,1);
% Loop through all taxonomic names to find AphiaID
for n_id = 1:num_taxa
  %fprintf('Finding AphiaID for %s\n',taxa.Name{n_id})
  % Query AphiaID 
  aphiaID = getAphiaID(obj,taxa.Name{n_id},args.marine_only);
  %   % See if parent name available
  %   if strcmp(taxa.Name_original,'Ctenophora(Ctenophora_XX)')
  %     keyboard
  %   end
  if flag_parent_available == 1
    % Make sure the entry is valid
    % for example, if using Ecotaxa output there are fields living with no
    % parent. 
    if ~isempty(taxa.Name_parent{n_id})
      % pull out parent AphiaID
      aphiaID_parent = getAphiaID(obj,taxa.Name_parent{n_id},args.marine_only);
      taxa.AphiaID_parent{n_id} = num2str(aphiaID_parent);
      % If query returns multiple matches, check previous aphiaID because
      % that may have already narrowed it down
      if isequal(aphiaID_parent,-999)
        if strcmp(taxa.Name{n_id-1},taxa.Name_parent{n_id}) && ~isequal(taxa.AphiaID{n_id-1},-999)
          taxa.AphiaID_parent{n_id} = taxa.AphiaID{n_id-1};
        end
      end
    else 
      taxa.Name_parent{n_id}    = '-';
      taxa.AphiaID_parent{n_id} = 'NULL';
    end
  end % IF flag_parent_available == 1
  
  %% Set AphiaID based on query return: NULL (none), or single match
  if isnan(aphiaID) % NO match (NULL)
    taxa.AphiaID{n_id} = 'NULL';
  elseif ~isequal(aphiaID,-999) % single match found 
    taxa.AphiaID{n_id} = num2str(aphiaID);
  end
  
  %% Set AphiaID based on query return: -999 (multiple) 
  % Special case because multiple matches were found. Need to go through
  % and verify choices. 
  if isequal(aphiaID,-999) % multiple matches found
    % Get one or more matching (max. 50) AphiaRecords for a given name.
    % Note: args.fuzzy is deprecated but included because required input
    aphiaIDs = getAphiaRecords(obj,taxa.Name{n_id},args.like, args.fuzzy, args.marine_only, args.offset);
    % ignore matches where the status is not accepted, e.g. "deleted"
    aphiaIDs = aphiaIDs(strcmp({aphiaIDs.status},'accepted'));
    % loop through matches until a match with parent, or grandparent is
    % found
    nmatch = 0; % Count matches
    flag_matchfound = 0;
    aphiaID_match   = '';
    while ~flag_matchfound
      % Update match index, do this at the beginning so only need to once
      % in case more scenarios are added later
      nmatch = nmatch + 1;
      % reached the end of the matches
      if nmatch > numel(aphiaIDs)
        % Set to NaN and continue to next taxa
        aphiaID_match.AphiaID  = NaN;
        aphiaID_match.lsid     = ''; % lsid : LifeScience Identifier. Persistent GUID for an AphiaID
        flag_matchfound = 1;
        continue
      end
      
      % See if scientific name is an exact (case insensitive) match
      % e.g. Cercozoa has 2 matches 'Cercozoa' and 'Cercozoa incertae sedis'
      if strcmpi(aphiaIDs(nmatch).scientificname,taxa.Name{n_id})
        %  Get the complete classification for one taxon. This also includes any sub or super ranks.
        aclass = getAphiaClassificationByID(obj,aphiaIDs(nmatch).AphiaID);
        classes  = struct();
        classes.ranks   = {};
        classes.sciname = {};
        % loop through all classes and pull out scientific names so can
        % verify that this will be the correct AphiaID. This is necessary
        % because of the format that getAphiaClassificationByID returns.
        while isstruct(aclass.child)
          classes.ranks   = [classes.ranks, aclass.rank];
          classes.sciname = [classes.sciname, aclass.scientificname];
          aclass = aclass.child;
        end
        % Check that parent is part of complete classification
        parent_idx = find(strcmp(taxa.Name_parent{n_id},taxa.Name));
        if any(strcmpi(taxa.Name_parent{n_id},classes.sciname))
          aphiaID_match   = aphiaIDs(nmatch);
          flag_matchfound = 1;
        elseif ~isempty(parent_idx)
          % Check that grandparent is part of complete classification
          % first, pull out index in table where parent name is the name
          % this is usually the one directly above...
          if any(strcmpi(taxa.Name_parent{parent_idx},classes.sciname))
            aphiaID_match   = aphiaIDs(nmatch);
            flag_matchfound = 1;
            %% No match for parent OR grandparent
          else
            % Check specific cases
            % Metazoa not in WoRMS... e.g. Ctenophora(Metazoa)
            % So, check if parent is Metazoa and next classificaiton above
            % is Animalia... in which case that is likely the match.
            if strcmp(taxa.Name_parent{n_id},'Metazoa') && strcmp(classes.sciname{end},'Animalia')
              taxa.AphiaID{n_id} = num2str(aphiaIDs(nmatch).AphiaID);
              aphiaID_match   = aphiaIDs(nmatch);
              flag_matchfound = 1;
            elseif strcmp(taxa.Name{n_id},'Ctenophora')
              if nmatch == numel(aphiaIDs)
                aphiaID_match.AphiaID  = NaN;
                aphiaID_match.lsid     = '';
                flag_matchfound        = 1;
              else
                %fprintf('%d/%d\n',nmatch,numel(aphiaIDs))
                continue
              end
            else % all other cases
              fprintf('Exact ScientificName match found for %s\n',taxa.Name{n_id})
              fprintf('.. BUT parent AND grandparent is not part of the taxa classificiaton\n')
              fprintf('.. Stop here to look at this case... entering keyboard mode\n')
              keyboard
            end
          end % IF grandparent is part of classification
        end % IF parent is part of classification
      end % IF scientific name matches the current taxa name
    end % LOOP through all matches until (hopefully) the correct AphiaID is found
    
    
    if flag_matchfound && isstruct(aphiaID_match)
      % Replace NaN with NULL
      if isnan(aphiaID_match.AphiaID)
        taxa.AphiaID{n_id} = 'NULL';
      else
        taxa.AphiaID{n_id} = num2str(aphiaID_match.AphiaID);
      end
      % Also pull out the Life Science Identifier (LSID), a persistent
      % globally unique identifier for biological objects in a
      % database, from WoRMS.
      taxa.scientificNameID{n_id} = aphiaID_match.lsid;
    else
      fprintf('Unexpected scenario... stopping here\n')
      keyboard
    end
  end
  % Replace NaN with NULL
  if strcmp(taxa.AphiaID_parent{n_id}, 'NaN')
    taxa.AphiaID_parent{n_id} = 'NULL';
  end
  % Set comment to empty string
  taxa.Comment{n_id} = '';
  
  % Build scientificName, when Name has a AphiaID match
  % An external, machine-readable and resolvable identifier, or object
  % number, that returns nomenclatural (not taxonomic) details of a name
  % (scientificNameID ) should also be included for each ROI.
  % A Life Science Identifier (LSID) used in a taxonomic database consists of
  % a uniform resource name (urn) that contains the following information in
  % order: network identifier, a root DNS name of the reference database, a
  % namespace, and the object number or taxon identifier (a living product)
  % that is unique to the biological object defined by that referenced
  % database.
  %  ---------------- Expected format -----------------------
  % [Network identifer]:[Reference database]:[Namespace]:[Taxon identifer]
  if ~strcmp(taxa.AphiaID{n_id}, 'NULL')
    % Check if ScientificNameID has already been assigned, do not overwrite
    if ~strncmp(taxa.scientificNameID{n_id},'urn',3)
      taxa.scientificNameID{n_id} = ['urn:lsid:marinespecies.org:taxname:' taxa.AphiaID{n_id}];
    end
  end
end %% LOOP THROUGH EACH ROW IN taxa TABLE

%% STEP4 | Check known special cases
% Add as necessary
if exist('special','var')
  fields = fieldnames(special);
  % Loop through known taxa names
  for nspec = 1:numel(fields)
    fname = fields{nspec};
    correct_name = fields{nspec};
    idx_correct  = strcmp(taxa.Name,correct_name);
    if all(idx_correct == 0)
      continue
    end
    % First, make sure it is a cell and not a char
    if ischar(special.(fields{nspec}))
      special.(fields{nspec}) = {special.(fields{nspec})};
    end
    % Loop through unknown taxa and match with the known name
    for nn = 1:numel(special.(fields{nspec}))
      special_name = special.(fields{nspec}){nn};
      idx_special  = strcmp(taxa.Name,special_name);
      
      % Only replace if it hasn't been identified yet
      if strcmp(taxa.AphiaID(idx_special),'NULL')
        if ~strcmp(taxa.AphiaID{idx_correct},'NULL')
          taxa.AphiaID(idx_special) = taxa.AphiaID(idx_correct);
        else
          aphiaID = getAphiaID(obj,correct_name,args.marine_only);
          % Only replace if NULL (no match), or -999 (multiple) were not returned
          if aphiaID > 0
            taxa.AphiaID{idx_special} = num2str(aphiaID);
            % Now check if special case is a parent anywhere, and if so, set that
            % AphiaID too.
            idx_parent_too = strcmp(taxa.Name_parent,taxa.Name{idx_special});
            if any(idx_parent_too)
              % Also make sure it hasn't already been ID - do not overwrite.
              if strcmp(taxa.AphiaID_parent(idx_parent_too),'NULL')
                taxa.AphiaID_parent{idx_parent_too} = taxa.AphiaID{idx_special};
              end
            end % Set AphiaID_parent for special cases
          end % Direct match found
        end % AphiaID has not been matched yet
      end % MATCH DIRECTLY, IF POSSIBLE
    end % LOOP through special KNOWN cases
  end % LOOP THROUGH TAXA NAMES THAT HAVE APHIAID TO MATCH SPECIAL CASES
end  % IF SPECIAL CASES EXIST

%% STEP5 | Match unknown taxa to their parent
% If parent is unknown, but still in the living categories, match to
% Eukaryota as suggested in OCB manual. 
idx_not_living = find(strcmp(taxa.Name,'not-living'));
if isempty(idx_not_living)
  idx_not_living = find(contains(taxa.Name_original,{'not-living' 'temporary'}));
end
for n_id = 1:num_taxa
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
      if numel(idx_not_living) > 1 && ismember(n_id,idx_not_living)
        continue
      elseif isempty(idx_not_living) || (~isempty(idx_not_living) && n_id < min(idx_not_living))
        taxa.scientificNameID{n_id} = 'urn:lsid:algaebase.org:taxname:86701';
        taxa.Comment{n_id} = 'miscellaneous living --> Eukaryota';
      end
    end % IF parent has, or does not have, matched AphiaID
  end % IF AphiaID has already been matched
end

%% STEP6 | MANUALLY Check unmatch cases by searching children of parent
if check_children_manual
  % For each given scientific name (may include authority), try to find one
  % or more AphiaRecords, using the TAXAMATCH fuzzy matching algorithm by
  % Tony Rees. This allows you to (fuzzy) match multiple names in one call.
  % Limited to 50 names at once for performance reasons
  for n_fuz = 1:num_taxa
    %  Only look at the unmatched cases where AphiaID_parent was matched
    if strcmp(taxa.AphiaID{n_fuz},'NULL') && ~strcmp(taxa.AphiaID_parent{n_fuz},'NULL') 
      parent_ID = str2double(taxa.AphiaID_parent{n_fuz});
      % Uses TAXAMATCH fuzzy matching algorithm by Tony Rees
      % Rees T (2014) Taxamatch, an Algorithm for Near (‘Fuzzy’) Matching
      % of Scientific Names in Taxonomic Databases. PLoS ONE 9(9): e107510.
      % https://doi.org/10.1371/journal.pone.0107510
      aphiaIDs = matchAphiaRecordsByNames(obj,taxa.Name{n_fuz},1);
      % If records are found, should return a structure
      if isstruct(aphiaIDs)
        fprintf('fuzzy matches found! stop here because need to set up scenario ...\n')
        keyboard
      elseif isfinite(parent_ID)
        % Get the direct children (max. 50) for a given AphiaID
        children = getAphiaChildrenByID(obj,parent_ID,1,1);
        % Only look at accepted entries
        children = children(strcmp({children.status},'accepted'));
        if isempty(children)
          continue
        end
        done_with_this_manual_check = 0; % allow user to stop script by entering 99
        while ~done_with_this_manual_check
          fprintf('\n')
          fprintf('No fuzzy matches found for %s\n',taxa.Name{n_fuz})
          fprintf('   %s(parent) --> %s(child)\n',taxa.Name_parent{n_fuz},taxa.Name{n_fuz})
          fprintf('   <0>  NONE (Default)\n')
          % Loop through all children
          for nchild = 1:numel(children)
            fprintf('   <%d>  %s\n',nchild,children(nchild).scientificname)
          end % loop through all children
          fprintf('            \n')
          fprintf('   <99> STOP\n')
          taxa_chc = input('   Enter correct number: ');
          if isempty(taxa_chc) || taxa_chc == 0
            fprintf('   ... okay, not selecting aphiaID\n')
            done_with_this_manual_check = 1;
          elseif taxa_chc == 99
            fprintf('   Entering keyboard mode, enter "dbcont" to continue\n')
            keyboard
          else
            taxa.AphiaID{n_fuz} = num2str(children(taxa_chc).AphiaID);
            taxa.Comment{n_fuz} = children(taxa_chc).scientificname;
            done_with_this_manual_check = 1;
          end
        end %% WHILE done_with_this_manual_check
      end
    end
  end % LOOP through taxa to try and match difficult cases
end % IF check_children_manual

%% STEP7 | Note all unmatched/non-conforming names in ScientificName [custom_namespace]
% The custom_namepsace is used to facilitate identification of
% non-conforming ROIs, custom definitions not found in a taxonomic
% authority must be provided in an external document file.
% For example, A Phytoplankton Taxonomy Working Group (“PTWG”) custom
% namespace was created to define several standardized names for common
% terms that are not currently defined by WoRMS or Algae Base. As of
% December 2020, this includes: ‘bad_image’, ‘bead’, ‘bubble’, ‘detritus’,
% ‘fecal_pellet’, and ‘other’. The term ‘other’ should only be used to
% describe a non-living particle.

% Loop through taxa table to find those note yet identified
for n_id = 1:num_taxa
  % First, check that scientificNameID has not already been assigned
  if strncmp(taxa.scientificNameID{n_id},'urn',3)
    continue
  end
  % Case 1: no match with CHILD or PARENT
  if strcmp(taxa.AphiaID{n_id},'NULL') && strcmp(taxa.AphiaID_parent{n_id},'NULL')
    taxa.scientificNameID{n_id} = ['namespace_unknown:' taxa.Name{n_id}];
    taxa.Comment{n_id} = 'non-conforming, search manually & add to namespace';
  % Case 2: no match with CHILD but match with PARENT, i.e. living
  elseif strcmp(taxa.AphiaID{n_id},'NULL') && ~strcmp(taxa.AphiaID_parent{n_id},'NULL')
    taxa.scientificNameID{n_id} = ['namespace_living:' taxa.Name{n_id}];
    taxa.Comment{n_id} = 'non-conforming, search manually & add to namespace';
  elseif ~strncmp(taxa.scientificNameID{n_id},'urn',3) % scientificName has not been set yet
    taxa.scientificNameID{n_id} = ['urn:lsid:marinespecies.org:taxname:' taxa.AphiaID{n_id}];
  end
end

%% STEP8 | Now that IDs have been identified where possible, pull out the AphiaName
for n_id = 1:num_taxa 
  if ~strcmp(taxa.AphiaID{n_id},'NULL')
    idnum = str2double(taxa.AphiaID{n_id});
    taxa.scientificName{n_id} = getAphiaNameByID(obj,idnum);    
  end
end

%% STEP9 | Search through all NULL living cases based on closest match in hierachy
if strcmp(data_provider,'ecotaxa')
  % ecotaxa will have "Name_original" field with family tree separated by >
  % For example....
  % Name_original=living>Eukaryota>Harosa>Rhizaria>Rhizaria X>RhizariaX sp.
  for n_id = 1:num_taxa
    if ismember(n_id,idx_not_living)
      continue
    end
    % If is living, but still has not been identified... search through
    % family tree and find match
    if strcmp(taxa.AphiaID{n_id},'NULL') && strcmp(taxa.AphiaID_parent{n_id},'NULL')
      family = strsplit(taxa.Name_original{n_id},'>');
      if numel(family) > 1
        family = fliplr(family); % search from lowest level
      end
      for nhierachy = 1:numel(family)
        family_name = family{nhierachy};
        if any(strcmp(taxa.scientificName,family_name))
          idx_match = find(strcmp(taxa.scientificName,family_name));
          if numel(idx_match) > 1
            idx_match = idx_match(1);
          end
          taxa.AphiaID(n_id) = taxa.AphiaID(idx_match);
          taxa.scientificName(n_id)   = taxa.scientificName(idx_match);
          taxa.scientificNameID(n_id) = taxa.scientificNameID(idx_match);
          taxa.Comment{n_id} = 'AphiaID --> matched within family hierarchy';    
          
        end
      end
    end
  end
end % ecotaxa

%% STEP11 | Save 
% save_filename = ['taxa_' date '.mat'];
% fprintf('Saving taxa matches to %s\n',[pwd filesep save_filename])
% save(save_filename,'taxa');

end %% MAIN FUNCTION








