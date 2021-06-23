function plotRandomRewardResponse(modelFile)
  
  %% Load data
  if ischar(modelFile)
    model             = load(modelFile);
  else
    model             = modelFile;
  end

  %% Find cells for which there were random reward trials
  cellData            = model.cellData(arrayfun(@(x) any([x.trials.condition_code] == 4), model.cellData));
  
  %%
  [pan,shape]         = makePanels(numel(cellData));
  
end
