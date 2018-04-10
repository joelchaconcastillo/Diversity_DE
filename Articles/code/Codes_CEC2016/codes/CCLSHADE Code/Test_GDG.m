function Test_GDG(fname,D)

lbvec = -100 * ones(1,D);
ubvec = 100 * ones(1,D);
epsilon = 1e-10;

Sep = [];
Non_Sep = [];
All_Groups = {};
FEVs = [];
num_groups = [];

for function_number = 1 : 30
    [seps, nonseps, FEs, DeltaMtx] = gdg(fname, function_number, D, lbvec, ubvec, epsilon);
    
    FEVs = [FEVs,FEs];
    
    Sep = [Sep,length(seps)>0];
    
    Non_Sep = [Non_Sep,numel(nonseps)>0];

    Groups = {};
    if(Sep(end)==0)
        for g = 1 : numel(nonseps)
            Groups{g} = nonseps{g};
        end
        num_groups(function_number) = numel(nonseps);
    elseif(numel(nonseps)==0)
        Groups = {seps};
        num_groups(function_number) = 1;
    else
        for g = 1 : numel(nonseps);
            Groups{g} = nonseps{g};
        end
        Groups = {seps Groups{1:end}};
        num_groups(function_number) = numel(nonseps) + 1;
    end
    
    All_Groups = {All_Groups{1:end} Groups};
    
    fprintf('\n\n Function %d Done\n\n',function_number);
end

ep = num2str(epsilon);

dimension = num2str(D);

filename = strcat('GDG_Data_',dimension,'_B','.mat');

save(filename,'Sep','Non_Sep','FEVs','num_groups','All_Groups');