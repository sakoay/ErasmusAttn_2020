function plotWaveletCoherence(continuousData)
  
  %% Define neuron and trials of interest
  neuron                = continuousData(strcmp({continuousData.cell_id},'190415.N3'));
  
  %% Construct behavioral event timing signal
  events                = struct();
  for what = fieldswithvalue(experiment.type, 'timing')'
    events.(what{:})    = accumfun(2, @(x) onehot(rmmissing(1 + x.(what{:})),x.duration,1,2), neuron.trials);
  end
  
  %%
  figure;   wcoherence(spikes, events.fix_on, milliseconds(1)); drawnow;
  figure;   wcoherence(spikes, events.saccade_onset, milliseconds(1)); drawnow;
  
end
