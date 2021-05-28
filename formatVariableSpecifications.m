function specs = formatVariableSpecifications(specs, conditions, abbreviate)
  
  %% Default arguments
  if ~exist('conditions', 'var')
    conditions  = [];
  end
  if ~exist('abbreviate', 'var') || isempty(abbreviate)
    abbreviate  = false;
  end

  %% Replace condition labels
  for iCond = 1:size(conditions,1)
    if abbreviate
      condStr   = regexprep(regexprep(conditions{iCond,2},' .*',''),'(rr).*','$1');
    else
      condStr   = conditions{iCond,2};
    end
    specs       = regexprep(specs, ['(condition_code=\S*)' conditions{iCond,1}], ['$1' condStr]);
  end
  
  %% Additional simplifications
  specs         = regexprep(specs, 'past_', 'past-');
  specs         = regexprep(specs, '_(direction|condition)[^=]*', '');
  specs         = regexprep(specs, '^(dir)[^=]*', '$1');
  specs         = regexprep(specs, '(_loc)ation[^=]*', '$1');
%   specs         = regexprep(specs, '^[^=]+=(.*[a-zA-Z].*)$', '$1');
  specs         = regexprep(specs, '^condition_code=', '');
  specs         = regexprep(specs, '_+', ' ');
  
end
