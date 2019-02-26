function CC_LSHADE_PopRestart_OneGroup(FUN,func_num,DIM,All_Groups,FEVs,num_groups,MaxFEVs)

lbound = -100;
ubound = 100;
lu = [lbound;ubound];

num_runs = 51;

T = [0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
%T = [0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
 
%mfes = [0.2,0.3,0.5];
%K = [0.2,0.5,1] * DIM;
mfes = [0.5,0.3,0.2];
K = [1,0.5,0.2] * DIM;
NumGroups = DIM./K;
max_nfes = MaxFEVs - FEVs(func_num);
    
CONV = [];

Regroup = 0;

for run = 1 : num_runs
    rand('state', sum(100*clock));
    T_Index = 1;
    conv = [];    

    % Create what every LSHADE suboptimizer needs
    % Fixed across all suboptimizers
    ipr = 18;%16
    p_best_rate = 0.11;%0.02;
    arc_rate = 1.4;%0.94;
    memory_size = 5;%7;
    
    % Randomly divide problem variables
    Groups = {};
    Kindex = 1;
    Problem_Variables = randperm(DIM);
    Num_Groups = NumGroups(Kindex);
    for group = 1 : Num_Groups
        Groups{end+1} = Problem_Variables(1+(group-1)*K(Kindex):group*K(Kindex));
    end
        
    % Different across suboptimizers
    problem_size = zeros(1,Num_Groups);
    for group = 1 : Num_Groups
        problem_size(group) = length(Groups{group});
    end
    pop_size = ipr * problem_size;
    max_pop_size = pop_size;
    min_pop_size = 4.0 * ones(1,Num_Groups);
    nfes = zeros(1,Num_Groups);
    memory_sf = 0.5 .* ones(memory_size, Num_Groups);%0.28
    memory_cr = 0.5 .* ones(memory_size, Num_Groups);%0.43
    memory_pos = ones(1,Num_Groups);
    %max_nfes_group = (max_nfes/Num_Groups) .* ones(1,Num_Groups);
    max_nfes_group = (mfes(Kindex)*max_nfes*pop_size(group)/sum(pop_size)) .* ones(1,Num_Groups);

    archive = struct('NP',round(arc_rate * pop_size(1)),'pop',zeros(0, problem_size(1)),'funvalues',zeros(0, 1));
    for group = 2 : Num_Groups
        archive = [archive;struct('NP',round(arc_rate * pop_size(group)),'pop',zeros(0, problem_size(group)),'funvalues',zeros(0, 1))];
    end
    
    % Initialize populations for all suboptimizers
    popold = {};
    pop = {};
    for group = 1 : Num_Groups
        popold{group} = lbound + (ubound-lbound) * rand(pop_size(group),problem_size(group));
        pop{group} = popold{group};
    end

    % Start with a randomly constructed context vector
    context_vector = zeros(1,DIM);
    for group = 1 : Num_Groups
        individual = ceil(rand*pop_size(group));
        context_vector(Groups{group}) = pop{group}(individual,:);
    end
    context_vector_fitness = feval(FUN,context_vector',func_num);
    max_nfes = max_nfes - 1;
    
    % Evaluate initial fitness for all populations of all suboptimizers 
    fitness = {};
    for group = 1 : Num_Groups
        population = repmat(context_vector,pop_size(group),1);
        population(:,Groups{group}) = pop{group};
        fitness{group} = feval(FUN, population', func_num);
        fitness{group} = fitness{group}';
        nfes(group) = nfes(group) + pop_size(group);
        
        [minfitness,minindex] = min(fitness{group});
        if(minfitness<context_vector_fitness)
            context_vector_fitness = minfitness;
            context_vector = population(minindex,:);
            %context_vector(Groups{group}) = pop{group}(minindex,:);
        end       
    end
    
    All_nfes = sum(nfes);
    
    if((All_nfes+FEVs(func_num))>=T(T_Index)*MaxFEVs)
        conv = [conv,context_vector_fitness];
        T_Index = T_Index + 1;
    end
    
    Stagnation = 0;
           
    % Cycling over LSHADE suboptimizers
    while(All_nfes<max_nfes)
        Improved = 0;
        
        % Time for a new group size
        if(All_nfes>(sum(mfes(1:Kindex))*max_nfes))
            % Randomly divide problem variables 
            % accoridng to the new group size
            Groups = {};
            Kindex = Kindex + 1;
            Problem_Variables = randperm(DIM);
            Num_Groups = NumGroups(Kindex);
            for group = 1 : Num_Groups
                Groups{end+1} = Problem_Variables(1+(group-1)*K(Kindex):group*K(Kindex));
            end           
                        
            % Different across suboptimizers
            problem_size = zeros(1,Num_Groups);
            for group = 1 : Num_Groups
                problem_size(group) = length(Groups{group});
            end
            pop_size = ipr * problem_size;
            max_pop_size = pop_size;
            min_pop_size = 4.0 * ones(1,Num_Groups);
            nfes = zeros(1,Num_Groups);
            memory_sf = 0.5 .* ones(memory_size, Num_Groups);
            memory_cr = 0.5 .* ones(memory_size, Num_Groups);
            memory_pos = ones(1,Num_Groups);
            %max_nfes_group = (max_nfes/Num_Groups) .* ones(1,Num_Groups);
            max_nfes_group = (mfes(Kindex)*max_nfes*pop_size(group)/sum(pop_size)) .* ones(1,Num_Groups);

            archive = struct('NP',round(arc_rate * pop_size(1)),'pop',zeros(0, problem_size(1)),'funvalues',zeros(0, 1));
            for group = 2 : Num_Groups
                archive = [archive;struct('NP',round(arc_rate * pop_size(group)),'pop',zeros(0, problem_size(group)),'funvalues',zeros(0, 1))];
            end

            % Initialize populations for all suboptimizers
            popold = {};
            pop = {};
            for group = 1 : Num_Groups
                popold{group} = lbound + (ubound-lbound) * rand(pop_size(group),problem_size(group));
                pop{group} = popold{group};
            end
            
            % Force current context vector into new populations
            for group = 1 : Num_Groups
                popold{group}(1,:) = context_vector(Groups{group});
                pop{group}(1,:) = popold{group}(1,:);
            end
            
            % Evaluate initial fitness for all 
            % populations of all suboptimizers 
            fitness = {};
            for group = 1 : Num_Groups
                population = repmat(context_vector,pop_size(group),1);
                population(:,Groups{group}) = pop{group};
                fitness{group} = feval(FUN, population', func_num);
                fitness{group} = fitness{group}';
                nfes(group) = nfes(group) + pop_size(group);
                All_nfes = All_nfes + pop_size(group);

                [minfitness,minindex] = min(fitness{group});
                if(minfitness<context_vector_fitness)
                    context_vector_fitness = minfitness;
                    context_vector = population(minindex,:);
                    %context_vector(Groups{group}) = pop{group}(minindex,:);
                end       
            end
            
            Stagnation = 0;
        end
        
        for group = 1 : Num_Groups
            pop{group} = popold{group};
            
            [temp_fit, sorted_index] = sort(fitness{group}, 'ascend');

            mem_rand_index = ceil(memory_size * rand(pop_size(group), 1));
            mu_sf = memory_sf(mem_rand_index,group);
            mu_cr = memory_cr(mem_rand_index,group);
            
            %% for generating crossover rate
            cr = normrnd(mu_cr, 0.1);
            term_pos = find(mu_cr == -1);
            cr(term_pos) = 0;
            cr = min(cr, 1);
            cr = max(cr, 0);
            
            %% for generating scaling factor
            sf = mu_sf + 0.1 * tan(pi * (rand(pop_size(group), 1) - 0.5));
            pos = find(sf <= 0);
            
            while ~ isempty(pos)
                sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
                pos = find(sf <= 0);
            end
            
            sf = min(sf, 1); 

            r0 = [1 : pop_size(group)];
            popAll = [pop{group}; archive(group).pop];
            [r1, r2] = gnR1R2(pop_size(group), size(popAll, 1), r0);
            
            pNP = max(round(p_best_rate * pop_size(group)), 2); %% choose at least two best solutions
            randindex = ceil(rand(1, pop_size(group)) .* pNP); %% select from [1, 2, 3, ..., pNP]
            randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
            pbest = pop{group}(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
            
            vi = pop{group} + sf(:, ones(1, problem_size(group))) .* (pbest - pop{group} + pop{group}(r1, :) - popAll(r2, :));
            vi = boundConstraintCC(vi, pop{group}, lu);
            
            mask = rand(pop_size(group), problem_size(group)) > cr(:, ones(1, problem_size(group))); % mask is used to indicate which elements of ui comes from the parent
            rows = (1 : pop_size(group))'; cols = floor(rand(pop_size(group), 1) * problem_size(group))+1; % choose one position where the element of ui doesn't come from the parent
            jrand = sub2ind([pop_size(group) problem_size(group)], rows, cols); mask(jrand) = false;
            ui = vi; ui(mask) = pop{group}(mask);
            
            children = repmat(context_vector,pop_size(group),1);
            children(:,Groups{group}) = ui;
            
            children_fitness = feval(FUN, children', func_num);
            children_fitness = children_fitness';
            nfes(group) = nfes(group) + pop_size(group);
            All_nfes = All_nfes + pop_size(group);
            
%             [minfitness,minindex] = min(children_fitness{group});
%             if(minfitness<context_vector_fitness)
%                 context_vector_fitness = minfitness;
%                 context_vector = children(minindex,:);
%                 %context_vector(Groups{group}) = ui(minindex,:);
%             end
               
            dif = abs(fitness{group} - children_fitness);
            
            %% I == 1: the parent is better; I == 2: the offspring is better
            I = (fitness{group} > children_fitness);
            goodCR = cr(I == 1);  
            goodF = sf(I == 1);
            dif_val = dif(I == 1);
            
            archive(group) = updateArchive(archive(group), popold{group}(I == 1, :), fitness{group}(I == 1));
            
            [fitness{group}, I] = min([fitness{group}, children_fitness], [], 2);
            
            popold{group} = pop{group};
            popold{group}(I == 2, :) = ui(I == 2, :);
            
            [minfitness,minindex] = min(fitness{group});
            if(minfitness<context_vector_fitness)
                context_vector_fitness = minfitness;
                context_vector(Groups{group}) = popold{group}(minindex,:);
                Improved = 1;
                Stagnation = 0;
            else
                Stagnation = Stagnation + pop_size(group);
            end

            num_success_params = numel(goodCR);
            
            if (num_success_params > 0) 
                sum_dif = sum(dif_val);
                dif_val = dif_val / sum_dif;

                %% for updating the memory of scaling factor 
                memory_sf(memory_pos(group),group) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);

                %% for updating the memory of crossover rate
                if ((max(goodCR) == 0) || (memory_cr(memory_pos(group),group))  == -1)
                  memory_cr(memory_pos(group),group)  = -1;
                else
                  memory_cr(memory_pos(group),group) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
                end

                memory_pos(group) = memory_pos(group) + 1;
                if (memory_pos(group) > memory_size)  
                    memory_pos(group) = 1; 
                end
            end
            
            %% for resizing the population size
            plan_pop_size = round((((min_pop_size(group) - max_pop_size(group)) / max_nfes_group(group)) * nfes(group)) + max_pop_size(group));

            if (pop_size(group) > plan_pop_size)
                reduction_ind_num = pop_size(group) - plan_pop_size;
                
                if (pop_size(group) - reduction_ind_num <  min_pop_size(group))
                    reduction_ind_num = pop_size(group) - min_pop_size(group);
                end

                pop_size(group) = pop_size(group) - reduction_ind_num;
                
                for r = 1 : reduction_ind_num
                  [valBest indBest] = sort(fitness{group}, 'ascend');
                  worst_ind = indBest(end);
                  popold{group}(worst_ind,:) = [];
                  pop{group}(worst_ind,:) = [];
                  fitness{group}(worst_ind,:) = [];
                end

                archive(group).NP = round(arc_rate * pop_size(group)); 

                if (size(archive(group).pop, 1) > archive(group).NP) 
                  rndpos = randperm(size(archive(group).pop, 1));
                  rndpos = rndpos(1 : archive(group).NP);
                  archive(group).pop = archive(group).pop(rndpos, :);
                end
            end
            
            if((All_nfes+FEVs(func_num))>=T(T_Index)*MaxFEVs)
                conv = [conv,context_vector_fitness];
                T_Index = T_Index + 1;
            end
    
            if(T_Index==length(T)+1)
                break;
            end
        end
        
        % Imrpovement after one evolutionary cycle
%         if(Improved==0)
%             Stagnation = Stagnation + 1;
%         else
%             Stagnation = 0;
%         end
        
        % Done when no improvement 
        % for 5 evolutionary cycles
        if(Stagnation>1000*DIM)
            if(Num_Groups==1)
                rgroup = 1;
                fmin = min(fitness{rgroup});
                fmax = max(fitness{rgroup});
                fconverged = ((fmax-fmin)<10^(-8));
                if(fconverged)
                    Stagnation = 0;
                    
                    %Different across suboptimizers
                    problem_size = zeros(1,Num_Groups);
                    problem_size(rgroup) = length(Groups{rgroup});
                    pop_size = ipr * problem_size;
                    max_pop_size = pop_size;
                    min_pop_size = 4.0 * ones(1,Num_Groups);
                    memory_sf = 0.5 .* ones(memory_size, Num_Groups);%0.28
                    memory_cr = 0.5 .* ones(memory_size, Num_Groups);%0.43
                    memory_pos = ones(1,Num_Groups);
                    %max_nfes_group = (max_nfes/Num_Groups) .* ones(1,Num_Groups);
                    max_nfes_group = (mfes(Kindex)*max_nfes*pop_size(rgroup)/sum(pop_size)-nfes(rgroup)) .* ones(1,Num_Groups);
                    nfes = zeros(1,Num_Groups);

                    archive = struct('NP',round(arc_rate * pop_size(1)),'pop',zeros(0, problem_size(1)),'funvalues',zeros(0, 1));

                    % Initialize populations for all suboptimizers
                    popold{rgroup} = lbound + (ubound-lbound) * rand(pop_size(rgroup),problem_size(rgroup));
                    pop{rgroup} = popold{group};

                    % Force current context vector into new populations
                    popold{rgroup}(1,:) = context_vector(Groups{rgroup});
                    pop{rgroup}(1,:) = popold{rgroup}(1,:);
                    
                    % Evaluate initial fitness for all 
                    % populations of all suboptimizers
                    fitness{rgroup} = {};
                    population = repmat(context_vector,pop_size(rgroup),1);
                    population(:,Groups{rgroup}) = pop{rgroup};
                    fitness{rgroup} = feval(FUN, population', func_num);
                    fitness{rgroup} = fitness{rgroup}';
                    nfes(rgroup) = nfes(rgroup) + pop_size(rgroup);
                    All_nfes = All_nfes + pop_size(rgroup);

                    [minfitness,minindex] = min(fitness{rgroup});
                    if(minfitness<context_vector_fitness)
                         context_vector_fitness = minfitness;
                         context_vector = population(minindex,:);
                         %context_vector(Groups{group}) = pop{group}(minindex,:);
                    end
                end
            else%if(Num_Groups>1)%>=0.1*mfes(Kindex)*max_nfes)
                %Stagnation = 0;

                % Collect parameters
                Population = zeros(pop_size(1),DIM);
                OldPopulation = zeros(pop_size(1),DIM);
                %Archive = zeros(archive(1).NP,DIM);
                for rgroup = 1 : Num_Groups
                    Population(:,Groups{rgroup}) = pop{rgroup};
                    OldPopulation(:,Groups{rgroup}) = popold{rgroup};
                    %Archive(:,Groups{group}) = archive(group).pop;
                end

                % Randomly redivide problem variables
                Groups = {};
                Problem_Variables = randperm(DIM);
                for rgroup = 1 : Num_Groups
                    Groups{end+1} = Problem_Variables(1+(rgroup-1)*K(Kindex):rgroup*K(Kindex));
                end

                % Redistribute parameters according
                % to the new decomposition
                for rgroup = 1 : Num_Groups
                    pop{rgroup} = Population(:,Groups{rgroup});
                    popold{rgroup} = OldPopulation(:,Groups{rgroup});
                    %archive(group).pop = Archive(:,Groups{group});
                end

                archive = struct('NP',round(arc_rate * pop_size(1)),'pop',zeros(0, problem_size(1)),'funvalues',zeros(0, 1));
                for rgroup = 2 : Num_Groups
                   archive = [archive;struct('NP',round(arc_rate * pop_size(rgroup)),'pop',zeros(0, problem_size(rgroup)),'funvalues',zeros(0, 1))];
                end

                Regroup = Regroup + 1;

                % Re-calculate fitness for all subpopulations
                for rgroup = 1 : Num_Groups
                    population = repmat(context_vector,pop_size(rgroup),1);
                    population(:,Groups{rgroup}) = popold{rgroup};
                    fitness{rgroup} = feval(FUN, population', func_num);
                    fitness{rgroup} = fitness{rgroup}';
                    nfes(rgroup) = nfes(rgroup) + pop_size(rgroup);
                    All_nfes = All_nfes + pop_size(rgroup);

                    [minfitness,minindex] = min(fitness{rgroup});
                    if(minfitness<context_vector_fitness)
                        context_vector_fitness = minfitness;
                        context_vector = population(minindex,:);
                        %context_vector(Groups{group}) = pop{group}(minindex,:);
                    end       
                end
            end
        end
    end
    
    % Saving the convergence 
    % results of the completed run
    if (run==1)
        CONV = conv;
    else
        if (length(conv)==size(CONV,2))
            CONV = [CONV;conv];
        elseif (length(conv)<size(CONV,2))
            conv = [conv,repmat(conv(end),1,size(CONV,2)-length(conv))];
            CONV = [CONV;conv];
        else
            CONV = [CONV,repmat(CONV(:,end),1,length(conv)-size(CONV,2))];
            CONV = [CONV;conv];
        end
    end

    fprintf('%d ',run);
    
%     % Randomly divide problem variables
%     % Done every run
%     Groups = {};
%     Problem_Variables = randperm(DIM);
%     for group = 1 : Num_Groups
%         Groups{end+1} = Problem_Variables(1+(group-1)*K:group*K);
%     end
end

CONV = CONV';

CONV = CONV - 100 * func_num;

number = num2str(func_num);
dimension = num2str(DIM);
filename = strcat('CCLSHADE_',number,'_',dimension,'.txt');
save(filename,'CONV','-ASCII');

fprintf('\n\n Function %d Done\t with Regrouping = %f times on average\n\n',func_num,Regroup/num_runs);