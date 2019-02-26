function CC_LSHADE_PopRestart(fid,FUN,func_num,DIM,All_Groups,FEVs,num_groups,MaxFEVs)

lbound = -100;
ubound = 100;
lu = [lbound;ubound];

num_runs = 5;

T = [0.01,0.02,0.03,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
%T = [0.0001,0.001,0.01,0.02,0.03,0.04,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0];
 

Groups = All_Groups{func_num};
max_nfes = MaxFEVs - FEVs(func_num);        
Num_Groups = num_groups(func_num);
%mfes = [0.2,0.2,0.2,0.2,0.2];
%mfes = [0.75,0.251];
    
CONV = [];

Restart = zeros(1,Num_Groups);
Width = 10^(-8) * (ubound-lbound);
    
for run = 1 : num_runs
    rand('state', sum(100*clock));
    T_Index = 1;
    conv = [];    
    index = 1;

    % Create what every LSHADE suboptimizer needs
    % Fixed across all suboptimizers
    p_best_rate = 0.11;
    arc_rate = 1.4;
    memory_size = 5;
    
    % Different across suboptimizers
    problem_size = zeros(1,Num_Groups);
    for group = 1 : Num_Groups
        problem_size(group) = length(Groups{group});
    end
    pop_size = 18 * problem_size;
    max_pop_size = pop_size;
    min_pop_size = 4.0 * ones(1,Num_Groups);
    nfes = zeros(1,Num_Groups);
    memory_sf = 0.5 .* ones(memory_size, Num_Groups);
    memory_cr = 0.5 .* ones(memory_size, Num_Groups);
    memory_pos = ones(1,Num_Groups);
    %max_nfes_group = (max_nfes/Num_Groups) .* ones(1,Num_Groups);
    max_nfes_group = (max_nfes*pop_size(group)/sum(pop_size)) .* ones(1,Num_Groups);
    %max_nfes_group = (mfes(index)*max_nfes*pop_size(group)/sum(pop_size)) .* ones(1,Num_Groups);

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
    % using the randomly constructed context vector
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
    
    %Stagnation = 0;
    Stagnation = zeros(1,Num_Groups);
           
    % Cycling over LSHADE suboptimizers
    while(All_nfes<max_nfes)
        Improved = 0;
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
                %Stagnation = 0;
                Stagnation(group) = 0;
            else
                %Stagnation = Stagnation + pop_size(group);
                Stagnation(group) = Stagnation(group) + pop_size(group);
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
            
            if(All_nfes>=T(T_Index)*max_nfes)
                conv = [conv,context_vector_fitness];
                T_Index = T_Index + 1;
            end
    
            if(T_Index==length(T)+1)
                break;
            end
            
            % Time to restart
            %if(All_nfes>(sum(mfes(1:index))*max_nfes))
            % No improvement for 1000D function evaluations
            % which is 10% of the total computational budget
            %if(Stagnation>=1000*DIM)
                %Stagnation = 0;
            % Check if any group has converged to restart it
            for rgroup = 1 : Num_Groups
                if(Stagnation(rgroup)>=1000*problem_size(rgroup))
                    %xmin = min(popold{rgroup});
                    %xmax = max(popold{rgroup});
                    %xdiff = find((xmax - xmin)<=Width);
                    %xconverged = (length(xdiff)==length(Groups{rgroup}));
                    fmin = min(fitness{rgroup});
                    fmax = max(fitness{rgroup});
                    fdiff = fmax - fmin;
                    fconverged = (fdiff<10^(-8));
                    % if all problem variables have converged
                    % or if all fitness are almost equal
                    if(fconverged)%&&xconverged)
                        %index = index + 1;
                        Restart(rgroup) = Restart(rgroup) + 1;
                        Stagnation(rgroup) = 0;
                                                
                        % Different across suboptimizers
                        pop_size(rgroup) = 18 * problem_size(rgroup);
                        max_pop_size(rgroup) = pop_size(rgroup);
                        min_pop_size(rgroup) = 4.0;
                        memory_sf(:,rgroup) = 0.5 .* ones(memory_size,1);
                        memory_cr(:,rgroup) = 0.5 .* ones(memory_size,1);
                        memory_pos(rgroup) = 1;
                        max_nfes_group(rgroup) = max_nfes_group(rgroup) - nfes(rgroup);
                        nfes(rgroup) = 0;
                        %max_nfes_group(rgroup) = ((max_nfes-All_nfes)*pop_size(rgroup)/sum(pop_size));
                        %max_nfes_group(group) = (mfes(index)*max_nfes*pop_size(group)/sum(pop_size)) .* ones(1,Num_Groups);

                        archive(rgroup) = struct('NP',round(arc_rate * pop_size(rgroup)),'pop',zeros(0, problem_size(rgroup)),'funvalues',zeros(0, 1));
                        
                        % Initialize populations for all suboptimizers
                        popold{rgroup} = lbound + (ubound-lbound) * rand(pop_size(rgroup),problem_size(rgroup));
                        pop{rgroup} = popold{rgroup};
                        
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
end

CONV = CONV';

CONV = CONV - 100 * func_num;

number = num2str(func_num);
dimension = num2str(DIM);
filename = strcat('CCLSHADE_',number,'_',dimension,'.txt');
save(filename,'CONV','-ASCII');

fprintf('\n\n Function %d Done with Restarting = ',func_num);
fprintf(fid,'Function %d \t',func_num);
for rgroup = 1 : Num_Groups
    fprintf('%0.2f ',Restart(rgroup)/num_runs);
    fprintf(fid,'%0.2f ',Restart(rgroup)/num_runs);
end
fprintf(fid,'\n');
fprintf('times on average\n\n');