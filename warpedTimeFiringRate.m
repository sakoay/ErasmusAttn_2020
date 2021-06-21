% WARPEDTIMEFIRINGRATE  Firing rate in warped time, which uses linear resampling between behavioral
% epochs to register all trials to the mean epoch duration across trials 
%
%  Inputs
% ========
%   trial         : 
%   eventName     : 
%   timeBin       : 
%
%  Outputs
% =========
%   firingRate    : 
%   sampleTime    : 
%
% Created :  20-Jun-2021 22:12:23
% Author  :  Sue Ann Koay (koay@princeton.edu)
%
function [firingRate, sampleTime] = warpedTimeFiringRate(trial, eventName, timeBin)

  %% Compute mean duration of intervals between events
  eventTimes              = accumfun(1, @(x) cellfun(@(y) x.(y), eventName(:)'), trial);
  eventTimes(:,end+1)     = [trial.duration];
  eventIntervals          = diff(eventTimes, 1, 2);
  assert( all(eventIntervals(:) > 0) );                                         % events must be occur in order from earliest to latest
  
  nominalTime             = [0, mean(eventIntervals, 1, 'omitnan')];
  
  %% Preallocate output
  nTimeBins               = ceil(sum(nominalTime) / timeBin);                   % N.B. last bin may have a smaller duration than timeBin; we will adjust for this later
  firingRate              = nan(numel(trial), nTimeBins);
  
  %% Linearly map spike times
  interp1
  
end
