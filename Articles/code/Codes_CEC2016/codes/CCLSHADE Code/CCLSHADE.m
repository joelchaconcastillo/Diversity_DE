function CCLSHADE(FUN)

for DIM = [10,30,50,100]        
    dimension = num2str(DIM);
    
    %filename = strcat('Diff_Grouping_Data_',dimension);
    
    filename = strcat('GDG_Data_',dimension);
    
    load(filename);
    
    filename = strcat('Regroup_Data_',dimension);
    
    fid = fopen(filename,'w');
    
    MaxFEVs = 200000;%10000 * DIM;

    tic; 
    for func_num = [18]
        % if only a small number of separable variables,
        % then add them to the smallest group
        if((Sep(func_num)==1)&&(numel(All_Groups{func_num}{1})<=0.1*DIM))            
            minimum = numel(All_Groups{func_num}{2});
            group = 2;
            for i = 3 : num_groups(func_num)
                if(numel(All_Groups{func_num}{i})<minimum)
                    minimum = numel(All_Groups{func_num}{i});
                    group = i;
                end
            end
            
            Sep_Variable = All_Groups{func_num}{1};
            All_Groups{func_num}{group} = [All_Groups{func_num}{group} Sep_Variable];
            All_Groups{func_num} = All_Groups{func_num}(2:end);
            num_groups(func_num) = num_groups(func_num) - 1;
        % if a large number of separable variables, 
        % then divide them into smaller groups
        elseif((Sep(func_num)==1)&&(numel(All_Groups{func_num}{1})>10))
            %All_Groups{func_num}
            %fprintf('%d\t%d\n',DIM,func_num);
            seps = All_Groups{func_num}{1};
            div_sepsize = 10;
            div_num = ceil(length(seps)/div_sepsize);
            div_seps = cell(1, div_num);
            div_startpts = 1:div_sepsize:length(seps);
            for i = 1:div_num-1
               div_seps{i} = seps(div_startpts(i):div_startpts(i+1)-1);
            end
            if (~isempty(seps))
               div_seps{div_num} = seps(div_startpts(div_num):end);
            end
            
           if((numel(div_seps{end})<=0.1*DIM)&&(numel(div_seps{end})<10))
               div_seps{end-1} = [div_seps{end-1} div_seps{end}];
               div_seps = div_seps(1:end-1);
               %div_seps
               %fprintf('GDG\t%d\n',func_num);   
           end
            
           All_Groups{func_num} = [All_Groups{func_num} div_seps];
           All_Groups{func_num} = All_Groups{func_num}(2:end);
           %All_Groups{func_num}
           num_groups(func_num) = num_groups(func_num) + numel(div_seps) - 1;
           %num_groups(func_num)
        end
        
        if(num_groups(func_num)>1)
            fprintf('GDG\t%d\n',func_num);
            CC_LSHADE_PopRestart(fid,FUN,func_num,DIM,All_Groups,FEVs,num_groups,MaxFEVs);
        else
            fprintf('Random\t%d\n',func_num);
            CC_LSHADE_PopRestart_OneGroup(FUN,func_num,DIM,All_Groups,FEVs,num_groups,MaxFEVs);
        end
    end
    toc/5
end
fclose(fid);
end
