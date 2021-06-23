% CONVERTTOGLMDATA    Convert Nico's data format to the one required by the Pillow lab neuroGLM
% code. This also performs some redefinitions and sanity checks.
%
%  Inputs
% ========
%   data    : array of data structures, with contents as in subset_cells of the Erasmus attention
%             project data from Nico.
%
%  Examples
% ==========
%   Plot data from a particular monkey:
%     [aniName, ~, aniIndex] = unique({subset_cells.monkey});
%     for iAni = 1:numel(aniName); convertToGLMData(subset_cells(aniIndex==iAni), aniName{iAni}, 'C:\Neuroscience\ErasmusAttn'); end
%
% Created :  05-Oct-2020 14:02:23
% Author  :  Sue Ann Koay (koay@princeton.edu)
%
function convertToGLMData(data, dataLabel, outputPath, timeBinMS)
  
  %% Default arguments
  if ~exist('timeBinMS', 'var') || isempty(timeBinMS)
    timeBinMS         = 50;
  end

  %% Setup experiment data structure 
  experiment          = buildGLM.initExperiment('ms', timeBinMS, dataLabel);
  
  % Ephys
  experiment          = buildGLM.registerSpikeTrain(experiment, 'SS', 'Simple spikes');
  
  % Task conditions
  experiment          = buildGLM.registerValue(experiment, 'trial_nr', 'Index of trial by temporal order in experiment');
  experiment          = buildGLM.registerValue(experiment, 'condition_code', '0 = aborted trial, 1 = error trial, 2 = correct trial, 3 = omitted reward, 4 = random reward');
  experiment          = buildGLM.registerValue(experiment, 'C_location', '1 = right, 3 = down, 5 = left, 7 = up');
  experiment          = buildGLM.registerValue(experiment, 'gap_direction', 'Orientation of the gap in the C, counterclockwise from north in degrees');
  experiment          = buildGLM.registerValue(experiment, 'past_condition', '0 = aborted trial, 1 = error trial, 2 = correct trial, 3 = omitted reward, 4 = random reward');
  experiment          = buildGLM.registerValue(experiment, 'past_C', '1 = right, 3 = down, 5 = left, 7 = top');
  experiment          = buildGLM.registerValue(experiment, 'past_gap', 'Orientation of the gap in the C');
  experiment          = buildGLM.registerValue(experiment, 'reward_duration', 'Amount of juice received');
  experiment          = buildGLM.registerValue(experiment, 'saccade_amplitude', 'Degrees of deflection in detected saccade');
  experiment          = buildGLM.registerValue(experiment, 'saccade_direction', 'Orientation of saccade, counterclockwise from north in degrees');
  experiment          = buildGLM.registerValue(experiment, 'past_saccade', 'Orientation of saccade, counterclockwise from north in degrees');

  % Task events
  experiment          = buildGLM.registerTiming(experiment, 'fix_on', 'Monkey initiates trial by fixation');
  experiment          = buildGLM.registerTiming(experiment, 'cue_start', 'Onset of priming cue for where on screen the C will appear');
  experiment          = buildGLM.registerTiming(experiment, 'cue_stop', 'Offset of priming cue');
  experiment          = buildGLM.registerTiming(experiment, 'c_start', 'Onset of the C that indicates the correct action direction');
  experiment          = buildGLM.registerTiming(experiment, 'c_stop', 'Offset of the C');
  experiment          = buildGLM.registerTiming(experiment, 'delay_start', 'Offset of the gray mask over the C');
  experiment          = buildGLM.registerTiming(experiment, 'choice_target', 'Onset of saccade target');
  experiment          = buildGLM.registerTiming(experiment, 'saccade_onset', 'Detected onset of saccade');
  experiment          = buildGLM.registerTiming(experiment, 'reward_start', 'Onset of juice delivery, or time at which trial was decided to be incorrect');
  
  % EyeLink tracking
  experiment          = buildGLM.registerContinuous(experiment, 'eye_xy', '(x,y) position of eye', 2);
%   experiment          = buildGLM.registerContinuous(experiment, 'pupil', 'Measure of pupil area', 1);     % FIXME : decide what to do with NaNs

  %% Loop over and register trials for each cell
  cellTrials          = cell(size(data));
  for iCell = 1:numel(data)
    %% UGLY : rename some fields to reduce confusion about definition
    data(iCell).trials= renamefield(data(iCell).trials, 'direction_C', 'C_location');
    data(iCell).trials= renamefield(data(iCell).trials, 'saccadeAmplitude', 'saccade_amplitude');
    
    %% Sanity checks
    if ~isfield(data(iCell).trials, 'n_minus_1')
      warning('convertToGLMData:trials', 'Data for %s cell %s (#%d) is missing n_minus_1 info, omitting this cell.', dataLabel, data(iCell).cell_id, iCell);
      continue;
    end
    
    invalidDuration   = find([data(iCell).trials.fixon2fixon] < 1);
    if ~isempty(invalidDuration)
      warning('convertToGLMData:trials', 'Data for %s cell %s (#%d) has trial(s) with zero duration: %s.', dataLabel, data(iCell).cell_id, iCell, mat2str(invalidDuration));
      data(iCell).trials(invalidDuration) = [];
    end
      
    %% HACK : skip some problematic data
    isBad             = ismember([data(iCell).trials.condition_code],1:2) & [data(iCell).trials.c_stop] < 0;
    if any(isBad)
      warning('convertToGLMData:c_stop', 'Invalid c_stop time for %s cell %s (#%d), trial(s) %s.', dataLabel, data(iCell).cell_id, iCell, mat2str(find(isBad)));
      data(iCell).trials(isBad) = [];
    end
    
    isBad             = any(cat(1,data(iCell).trials.event_presentations) > 1, 2);
    if any(isBad)
      warning('convertToGLMData:event_presentations', 'Multiple presentations of events for %s cell %s (#%d), trial(s) %s.', dataLabel, data(iCell).cell_id, iCell, mat2str(find(isBad)));
      data(iCell).trials(isBad) = [];
    end
    
    isBad             = any(cat(1,data(iCell).trials.event_presentations) > 1, 2);
    if any(isBad)
      warning('convertToGLMData:event_presentations', 'Multiple presentations of events for %s cell %s (#%d), trial(s) %s.', dataLabel, data(iCell).cell_id, iCell, mat2str(find(isBad)));
      data(iCell).trials(isBad) = [];
    end
    
    isBad             = arrayfun(@(x) isempty(x.SS), data(iCell).trials);
    if any(isBad)
      warning('convertToGLMData:SS', 'No spikes for %s cell %s (#%d), trial(s) %s.', dataLabel, data(iCell).cell_id, iCell, mat2str(find(isBad)));
      data(iCell).trials(isBad) = [];
    end
    
    %% Convert trial data to GLM format
    cellTrials{iCell} = mstruct(size(data(iCell).trials), [{'duration'}; fieldnames(experiment.type)]);
    for iTrial = 1:numel(data(iCell).trials)
      %% Define source data and any extra computed fields here
      source          = data(iCell).trials(iTrial);
      source.n_minus_1= renamefield(source.n_minus_1, 'direction_C', 'C_location');
      
      source.delay_start      = source.c_stop + 100;
      if isempty(source.saccade_landmarks)
        source.saccade_onset  = nan;
      else
        source.saccade_onset  = source.saccade_landmarks(1);
      end
      
      %% HACK : reinterpret trial start as "fixation" start for random reward trials
      if source.condition_code == 4
        %% Adjust timings of behavioral events
        for what = fieldswithvalue(experiment.type, 'timing')'
          if source.(what{:}) >= 0
            source.(what{:})  = double(source.(what{:})) - double(source.trial_start);
          end
        end
        source.fixon2fixon    = double(source.fixon2fixon) - double(source.trial_start);
        source.fix_on         = 0;
        
        %% Adjust spike timings
        for what = fieldswithvalue(experiment.type, 'spike train')'
          source.(what{:})    = source.(what{:}) - source.trial_start/1000;   % in seconds
        end
      end
      
      %% Flatten past-trial fields for convenience
      for what = fieldnames(source.n_minus_1)'
        source.(['past_', regexprep(what{:},'_.*','')]) = source.n_minus_1.(what{:});
      end
      
      %% Deduce saccade direction
      switch source.condition_code
        case 1        % error trial
          source.saccade_direction  = mod(source.gap_direction + 180, 360);
        case 2        % correct trial
          source.saccade_direction  = source.gap_direction;
        otherwise
          source.saccade_direction  = nan;
      end
      
      switch source.past_condition
        case 1        % error trial
          source.past_saccade       = mod(source.past_gap + 180, 360);
        case 2        % correct trial
          source.past_saccade       = source.past_gap;
        otherwise
          source.past_saccade       = nan;
      end
      
      %% Truncate data to only the extent of this trial
      trialDuration   = double(source.fixon2fixon);
      trialRange      = 4000 + [0, trialDuration - 1];
      if trialRange(2) > size(source.eyeX,1)
        %% FIXME : pad with NaNs for now
        warning('convertToGLMData:trialRange', 'Range of Ephys recordings (%.4g ms) exceeds available behavioral data, for %s cell %s (#%d), trial %d.', trialDuration, dataLabel, data(iCell).cell_id, iCell, iTrial);
        range         = trialRange(1):size(source.eyeX,1);
        source.eye_xy = [source.eyeX(range), source.eyeY(range)];
        source.eye_xy(end+1:trialRange(2) - trialRange(1) + 1,:)  = nan;
      else
        range         = trialRange(1):trialRange(2);
        source.eye_xy = [source.eyeX(range), source.eyeY(range)];
      end
      

      %% Add all prespecified data of interest
      trial           = buildGLM.newTrial(experiment, trialDuration);
      info            = fieldnames(experiment.type);
      for iField = 1:numel(info)
        switch experiment.type.(info{iField})
          case 'spike train'
            %% Convert spike timings from seconds to ms, and select only spikes within the range of this trial
            ephys     = double(source.(info{iField})) * 1000;       % convert s -> ms
            ephys( ephys < 0 | round(ephys) >= trialDuration ) = [];

            %% Align timing info to trial start instead of 0 = fixation on
            trial.(info{iField})     = ephys;
          case 'timing'
            %% Align timing info to trial start instead of 0 = fixation on
            if source.condition_code == 4 && source.(info{iField}) < 0
              %% Special case for non-trials
              trial.(info{iField})   = nan;
            else
              trial.(info{iField})   = double(source.(info{iField}));
            end
          case 'continuous'
            %% Compute temporally binned average
            trial.(info{iField})     = rebin(double(source.(info{iField})), experiment.binSize, 1);
          otherwise
            %% Verbatim
            trial.(info{iField})     = double(source.(info{iField}));
        end
      end
      
      %% Check validity of trial structure, and store separately for this cell
      buildGLM.addTrial(experiment, trial, 1);
      cellTrials{iCell}(iTrial) = rmfield(trial, 'expt');       % reduce output size
    end
  end
  
  %% Retain some identification info for each cell
  cellData            = keepfield(data, '^(monkey|cell_id)$');
  [cellData.trials]   = cellTrials{:};
  cellData(cellfun(@isempty, cellTrials)) = [];
  
  %% Determine output location
  outputFile          = fullfile(outputPath, sprintf('glmExpt_%s_%dcells_%dms.mat', dataLabel, numel(cellData), timeBinMS));
  
  %% Save output
  fprintf(' ======>  %s\n', outputFile);
  save(outputFile, 'experiment', 'cellData');

end
