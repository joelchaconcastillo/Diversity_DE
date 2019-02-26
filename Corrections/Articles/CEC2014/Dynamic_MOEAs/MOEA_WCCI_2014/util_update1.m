function util_update1()
%Search Utility updation.
%   This function update the subproblem's search utility, using the
%   improvement of the subproblems' objective vaue since the last caculation.
    
    global subproblems idealpoint;
    shrehold= 0.01;
    
    %allweight = [subproblems.weight];
    %curps = [subproblems.curpoint];
    %oldps = [subproblems.oldpoint];
    nsize=length(subproblems);
    util=ones(nsize,1);
    newobj=[];oldobj=[];
  objdim=length(subproblems(1).curpoint.objective);
  for i=1:nsize
      
      oldobj=[oldobj; subproblems(i).oldpoint.objective'];
  end
  for i=1:nsize
      
      newobj=[newobj; subproblems(i).curpoint.objective'];
  end
  utilmat=zeros(nsize,objdim);
  %maxobj=max(nobj);
  %minobj=min(nobj);
  
  %p=(maxobj-minobj)/10 ;    
  
  
    
    for  i=1:nsize
      
      
      for j=1:objdim
          
          if (newobj(j)-oldobj(i,j)<=-shrehold)
              utilmat(i,j)=1;
          else
             utilmat(i,j)=0.99-(newobj(j)-oldobj(i,j));
          end
      end
      
    end
    util_prod=prod(utilmat,2);
    
    util=util_prod;%util.*util_prod;
    for i=1:nsize
        subproblems(i).utility=util(i);
    end
    [subproblems.oldpoint] = subproblems.curpoint;
    %newobj=subobjective(allweight, [curps.objective], idealpoint, 'te');
    %oldobj=subobjective(allweight, [oldps.objective], idealpoint, 'te' );
    
    %delta =  oldobj - newobj;
    
    %in rare case the delta can be less then 0. change it in that
    %situition.
    %[subproblems(delta<0).curpoint] = subproblems(delta<0).oldpoint;
    %[subproblems(delta<0).utility] = deal(1);
    
    %[subproblems(delta>=shrehold).utility] = deal(1.0);
    
    %update = and(delta<shrehold, delta>=0);
    %update = and(delta<shrehold, delta>=0);
    %update = delta<shrehold;
    % make the updation.
    %if any(update)
     %   util = (0.95 + 0.05*delta(update)/shrehold) ... 
      %      .*[subproblems(update).utility];
        %util = min(util, 1.0);
       % cellutil = num2cell(util);
        %[subproblems(update).utility] = cellutil{:};
    %end
    
    %back up the old optimal values.
    
end
