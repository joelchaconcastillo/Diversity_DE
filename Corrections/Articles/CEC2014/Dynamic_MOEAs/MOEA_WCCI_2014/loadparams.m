function params = loadparams(mop, propertyArgIn)
%LOADPARAMS Load the parameter settings from external file.
% file format 
%   Detailed explanation goes here

    params = get_structure('parameter');
    % set the default!!

    % General Setting for MOEA/D
    params.seed = 177; % random seed.
    params.popsize = 300; % population size.
    params.niche = 20;  % neighbourhhood size
    params.dmethod = 'ts'; % decomposition method of choice.
    params.iteration = 300-1; % total iteration.=window x 60
    params.evaluation = 30000; % total function evaluation number=iteration x particle.

    % new MOEA/D setting.
    params.updateprob = 0.9; % percantege to update the neighbour or the whole population.
    params.updatenb = 2; % the maximum updation number.

    % DE setting.
    params.F = 0.5; % the F rate in DE.
    params.CR = 1; % the CR rate in DE.

    % R-MOEA/D specific setting.
    params.dynamic = 0; % weather to use R-MOEA/D or a normal MOEA/D, 1 for R-MOEA/D
    params.selportion = 5; % the selection percentage of in R-MOEA/D
    params.decayrate = 0.95; % the decay rate in R-MOEA/D

    % set for the default for a specific problem here.
    switch lower(mop.name)
        case {'fda1','zjz','zdt11','fda2','zjz2','fda3'}
            params.popsize = 100;
            params.niche = 10;
            %params.updatenb=6;
        case {'cf1_var1','cf1_var2','dyn1','dyn1_ps','dyn1_pf','ud1','ud2','ud3','ud4','ud5','ud7','ud8','ud9','ud10','ud11','ud13','ud14'}
            params.popsize = 300;
            params.niche = 30;
        case {'ud6','ud12'}
            params.popsize = 500;
            params.niche = 50;
        otherwise
    end

    % handle the parameters passed in from the function directly!
    % it has priority higher than the default values.
    while length(propertyArgIn)>=2
        prop = propertyArgIn{1};
        val=propertyArgIn{2};
        propertyArgIn=propertyArgIn(3:end);

        switch prop
            case 'popsize'
                params.popsize=val;
            case 'niche'
                params.niche=val;
            case 'evaluation'
                params.evaluation=val;
            case 'selportion'
                params.selportion=val;                
            case 'dynamic'
                params.dynamic=val;
            case 'seed'
                params.seed=val;                
            otherwise
                warning('moea doesnot support the given parameters name');
        end
    end

end

