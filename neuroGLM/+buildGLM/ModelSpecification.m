classdef ModelSpecification < handle
  
  %_________________________________________________________________________________________________
  properties (SetAccess = public)
    expt        = []
    covar       = []
    idxmap      = []
    edim        = []
    stimLabel   = []
  end
  
  %_________________________________________________________________________________________________
  methods
    
    %----- Constructor 
    function obj = ModelSpecification(expt)
      if nargin > 0
        obj.expt  = expt;
      end
    end
    
  end
  
end
