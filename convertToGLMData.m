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

  %% Determine output location
  outputFile          = fullfile(outputPath, sprintf('glmExpt_%s_%dcells_%dms.mat', dataLabel, numel(data), timeBinMS));
  
  %% Setup experiment data structure 
  experiment          = buildGLM.initExperiment('ms', timeBinMS, dataLabel);
  
  % Ephys
  experiment          = buildGLM.registerSpikeTrain(experiment, 'SS', 'Simple spikes');
  
  % Task conditions
  experiment          = buildGLM.registerValue(experiment, 'trial_nr', 'Index of trial by temporal order in experiment');
  experiment          = buildGLM.registerValue(experiment, 'condition_code', '0 = aborted trial, 1 = incorrect trial, 2 = correct trial, 3 = omitted reward, 4 = random reward');
  experiment          = buildGLM.registerValue(experiment, 'direction_C', '1 = right, 3 = down, 5 = left, 7 = top');
  experiment          = buildGLM.registerValue(experiment, 'gap_direction', 'Orientation of the gap in the C');
  experiment          = buildGLM.registerValue(experiment, 'reward_duration', 'Amount of juice received');
  experiment          = buildGLM.registerValue(experiment, 'saccade_amplitude', 'Degrees of deflection in detected saccade');

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
  parfor iCell = 1:numel(data)
    cellTrials{iCell} = mstruct(size(data(iCell).trials), [{'duration'}; fieldnames(experiment.type)]);
    
    for iTrial = 1:numel(data(iCell).trials)
      %% Define source data and any extra computed fields here
      source          = data(iCell).trials(iTrial);
      source.delay_start      = source.c_stop + 100;
      if isempty(source.saccade_landmarks)
        source.saccade_onset  = nan;
      else
        source.saccade_onset  = source.saccade_landmarks(1);
      end
      
      %% Truncate data to only the extent of this trial
      % We should restrict the data to only the Ephys recording period because we have no info about
      % neural activity outside of this window, and should not assume that it is zero.
      % We have the trial start time according to the EyeLink clock, and need to convert it to the
      % behavioral clock which is per-trial. To do this we use the reward time for which we have
      % both EyeLink and behavioral timestamps.
      
      % Start of trial relative to reward time
      startRelReward  = double(source.timestamp_el_reward) - double(source.ephys_start);
      % Start of Ephys recording according to the behavioral clock 
      ephysStart      = double(source.reward_start) - double(startRelReward);
      % End of Ephys recording according to the behavioral clock 
      ephysStop       = ephysStart + double(source.triallength) - 1;
      
      %% Apply truncation window to 1ms sampled behavioral data 
      trialDuration   = ephysStop - ephysStart + 1;
      trialRange      = 4000 + ephysStart + [0, trialDuration - 1];
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
            ephys( ephys < ephysStart | ephys > ephysStop ) = [];

            %% Align timing info to trial start instead of 0 = fixation on
            trial.(info{iField})     = ephys - ephysStart;
          case 'timing'
            %% Align timing info to trial start instead of 0 = fixation on
            if source.condition_code == 4
              %% Special case for non-trials
              trial.(info{iField})   = nan;
            else
              trial.(info{iField})   = double(source.(info{iField})) - ephysStart;
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
  
  %% Save output
  fprintf(' ======>  %s\n', outputFile);
  save(outputFile, 'experiment', 'cellData');

end
