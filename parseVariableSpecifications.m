function specs = parseVariableSpecifications(specs, lookup, abbreviate)
  %% Default arguments
  if ~exist('abbreviate', 'var') || isempty(abbreviate)
    abbreviate  = false;
  end
  
  %% Parse description string for specifications
  specs         = regexp(specs, '([0-9]+)\s*=\s*([^,]+)', 'tokens');
  specs         = catcell(1, specs, 2);
  
  %% Return labels for the given set of values
  if nargin < 2 || isempty(lookup)
    if ~isempty(specs)
      specs     = specs{:};
    end
  elseif ~isempty(specs)
    %% Additional simplifications
    specs       = specs{:};
    if abbreviate
      specs(:,2)= regexprep(specs(:,2), '^(.)(\S+)$', '${upper($1)}');
      specs(:,2)= regexprep(specs(:,2), '\s.*', '');
    end
  
    %% Convert indices to labels
    values      = cellfun(@str2num, specs(:,1));
    if iscell(lookup)
      lookup    = cellfun(@str2num, lookup, 'UniformOutput', false);
      [~,iMatch]= cellfun(@(x) ismember(x,values), lookup, 'UniformOutput', false);
      specs     = cellfun(@(x) ['[' strjoin(specs(x,2),';') ']'], iMatch, 'UniformOutput', false);
      return;
    elseif ischar(lookup)
      [~,iMatch]= ismember(str2num(lookup), values);
    else
      [~,iMatch]= ismember(lookup, values);
    end
    specs       = specs(iMatch,2);
  elseif iscell(lookup) || ischar(lookup)
    specs       = lookup;
  else
    specs       = arrayfun(@num2str, lookup, 'uniformoutput', false);
  end
end
