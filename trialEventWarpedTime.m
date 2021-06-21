% TRIALEVENTWARPEDTIME  Linear resampling between behavioral epochs to register all trials to the
% mean epoch duration across trials 
%
%  Inputs
% ========
%   experiment    : GLM experiment
%   eventName     : Ordered list of events between which to linearly warp time
%
%  Outputs
% =========
%   origTime      : Original time coordinate for all trials, relative to eventName{1}
%   warpedTime    : Warped time coordinate for all trials
%
% Created :  20-Jun-2021 16:35:48
% Author  :  Sue Ann Koay (koay@princeton.edu)
%
function [origTime, warpedTime] = trialEventWarpedTime(experiment, eventName)

  %% Compute mean duration of intervals between events
  eventTimes              = accumfun(1, @(x) cellfun(@(y) experiment.binfun(x.(y)), eventName(:)'), experiment.trial);
  eventTimes(:,end+1)     = cellfun(@(x) numel(x)+1, experiment.trialBins);     % add an event that starts after end of trial, for trial duration
  eventIntervals          = diff(eventTimes, 1, 2);
  assert( all(eventIntervals(:) > 0) );                                         % events must be occur in order from earliest to latest
  
  nominalTime             = [0, mean(eventIntervals, 1, 'omitnan')];
  
  %% Define the original time coordinate starting from the first event for each trial
  origTime                = accumfun(2, @(x) x - x(1), experiment.trialBins);

  %% Define warped time with linear interpolation between events
  warpedTime              = accumfun(2, @(x) accumfun(2, @(y) butlast(linspace(nominalTime(y), nominalTime(y+1), eventIntervals(x,y)+1)), 1:size(eventIntervals,2)), 1:size(eventTimes,1));
  
end
