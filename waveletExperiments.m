function waveletExperiments(cellData, experiment)
  
  %% Define neuron and trials of interest
  neuron                = cellData(strcmp({cellData.cell_id},'190415.N3'));
  selTrials             = ismember([neuron.trials.condition_code], 2);
  [sum(selTrials), mean(selTrials)]
  neuron.trials         = neuron.trials(selTrials);
  neuron.trials(2:2:end)= [];
  spikes                = accumfun(2, @(x) onehot(1+round(x.SS),x.duration,1,2), neuron.trials);
  
  %% Construct behavioral event timing signal
  events                = struct();
  for what = fieldswithvalue(experiment.type, 'timing')'
    events.(what{:})    = accumfun(2, @(x) onehot(rmmissing(1 + x.(what{:})),x.duration,1,2), neuron.trials);
  end
  
  %%
  figure;   wcoherence(spikes, events.fix_on, milliseconds(1)); drawnow;
  figure;   wcoherence(spikes, events.saccade_onset, milliseconds(1)); drawnow;
  
end
