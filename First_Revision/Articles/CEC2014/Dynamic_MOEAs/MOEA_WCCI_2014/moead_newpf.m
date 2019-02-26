function [pareto metric1 p1 p2 p3 p4 p5 ]= moead_newpf( mop, varargin)  %done with PF forward march with bound checking
%MOEAD run moea/d algorithms for the given mop.
% MOP could be obtained by function 'testmop.m'.
% the controlling parameter can be passed in by varargin, the folliwing
% parameters are defined in here. More other parameters can be passed by
% modify loadparams.m problem by problem.
%   seed: the random seed.
%   popsize: The subproblem's size.
%   niche: the neighboursize, must less then the popsize.
%   evaluation: the total evaluation of the moead algorithms before finish.
%   dynamic: whether to use dynamic resource allocation.
%   selportion: the selection portion for the dynamic resource allocation

    %global variable definition.
    global subproblems params itrCounter evalCounter rnduni window1 step1 fnum at bt ct dt et Fps Gps Hpf Mpf Gpf;
    %global idealpoint objDim parDim evalCounter;
    
    %load the parameters.
    params=loadparams(mop, varargin);
    
    %set the random generator.
    %rs = RandStream.create('mt19937ar', 'Seed', params.seed);
    %RandStream.setDefaultStream(rs);
    %seed = rem((params.seed+23),1377);
    %rnduni = -seed;
   
    %the counters.
    evalCounter = 0;
    itrCounter = 1;
    
    %and Initialize the algorithm.
    init(mop,itrCounter);
    pareto=[subproblems.curpoint];
    archive_present=[pareto.objective];
    archive_past=[pareto.objective];
    parameter_present=[pareto.parameter];
    parameter_past=[pareto.parameter];
    metric1=[];
    
 while ~terminate()
        
            dim=mop.pd;
         num=1;
          if mod(itrCounter,window1)==0
         r=randperm(5);
         fnum=r(1);
    switch fnum
        case {1}
        G1=sin(0.4*pi*at);
        F=ceil(dim*G1);
        Fps=F*ones(dim-1,num);
        at=at+1;       
        case {2}        
        G2=sin(0.4*pi*bt);
        Gps=G2*ones(1,num);
        bt=bt+1;
        case {3}    
            Hpf=0.5+abs(sin(0.4*pi*ct));
            ct=ct+1;
        case {4}
            Mpf=0.5+abs(sin(0.4*pi*dt));
            dt=dt+1;
        case {5}
            Gpf=abs(sin(0.4*pi*et));
            et=et+1;
    end
    
    end

     
     
        evolve(mop,itrCounter); % one generation of evaluation.
        
        pareto=[subproblems.curpoint];
        front=[pareto.objective];
        if itrCounter==229
            p1=front;
        end
        if itrCounter==234
            p2=front;
        end
        
        if itrCounter==239
            p3=front;
        end
        
        if itrCounter==244
            p4=front;
        end
        
        if itrCounter==249
            p5=front;
        end
        
        
        set=[pareto.parameter];
         %metric=igd2(front, mop.name, mop.pd,itrCounter) ;
        metric=igd2(front, mop.name, mop.pd,itrCounter) ;
        % metric=efmetric(set, mop.name, mop.pd,itrCounter) ;
         metric1=[metric1 metric];
         disp(sprintf('Itr:%d\tMetric:%1.6f', ...
        itrCounter, metric));
    itrCounter=itrCounter+1;
%         if (rem(itrCounter,50)==0) % updating of the utility.
%             util_update(); 
%             status();
%         end  
%% basic operator for forward march


           
        if mod(itrCounter,window1)==0
           dd=[];
            archive_past1=archive_past;
            archive_past=archive_present;
            archive_present=front;
            parameter_past1=parameter_past;
            parameter_past=parameter_present;
            parameter_present=[pareto.parameter];
            blank=[];
            no=length(archive_past);
            for i=1:no
                p=archive_present(:,i);
                tempdist=[];
                for j=1:no
                    q = archive_past(:,j);
                    r=(p-q).^2;
                    s=sum(r);
                    s=sqrt(s);
                    tempdist=[tempdist s];
                end 
                [Y, M]=min(tempdist);
                
                p=archive_past(:,M);
                tempdist1=[];
                for j=1:no
                    q = archive_past1(:,j);
                    r=(p-q).^2;
                    s=sum(r);
                    s=sqrt(s);
                    tempdist1=[tempdist1 s];
                end 
                 D= tempdist1==0;
                 tempdist1(D)=inf;
              
                [Y2, M1]=min(tempdist1);
                
                 P1=parameter_present(:,i)-parameter_past(:,M);
                 P2=parameter_past(:,M)-parameter_past1(:,M1);
                 k=abs(P1)-abs(P2);
                 
                  for j=1:mop.pd
%                    if P1(j)<=P2(j)
%                       epsilon(j)=-k(j)/1;
%                    else
%                      epsilon(j)=k(j)/1;
%                    end
                        x=k(j);
                        
%                      epsilon(j)=1+tanh(x);%+abs(randn);
                       if rand <0.5
                           epsilon(j)=1+tanh(x);%+abs(randn);
                        else
                         epsilon(j)=1+(x)*abs(randn);
                        end
                  end
                  m=P1;
                 % dd=[dd;k(j) epsilon(j)-1]
%                  m=mean(P1+P2,2);
%                   sd=std([P1 P2],0,2)/40;
%                  ltemp=length(m);
                  pipi=m.*epsilon';%sd.*randn(ltemp,1);
                  %newposition=parameter_present(:,i)+(parameter_present(:,i)-parameter_past(:,M));   %% change has to be done
                  %newposition=(3*(parameter_present(:,i)-parameter_past(:,M))+parameter_past1(:,M1));
                  newposition=parameter_present(:,i)+pipi;%(P1+P2)/2;
                  bound=mop.domain;
                  c=0;
                  for j=1:mop.pd
                      if (newposition(j))<bound(j,1)||(newposition(j))>bound(j,2)
                          newposition(j)=parameter_present(j,i);%bound(j,1)+(bound(j,2)-bound(j,1))*rand;
                          c=c+1;
                      end
                  end
%                   blank=[blank c];
                    sigma=(Y^2)/(4*mop.pd);
                    epsilon1=sigma*randn;
%                     dist=(parameter_present(:,i)-parameter_past(:,M)).^2;
%                     epsilon=(sum(dist))^0.5;
%                 sigma=1.5*epsilon*randn(mop.pd,1);
%                 subproblems(i).curpoint.parameter=parameter_present(:,i)+sigma;
                   subproblems(i).curpoint.parameter=newposition;%+epsilon1*ones(mop.pd,1);
                
                ind=evaluate(mop,subproblems(i).curpoint,itrCounter);
                subproblems(i).curpoint.objective=ind;
            end
%       length(find(blank))
%      blank
        end
 end
    %display the result.subproblems(i).optimal=newobj(i);
    
    %pp=[pareto.objective];
    %scatter(pp(1,:), pp(2,:));
    %disp(sprintf('total time used %u', etime(clock, starttime)));
 end

% The evoluation setp in MOEA/D
function evolve(mop,itrCounter)
  global subproblems idealpoint params rnduni window1 step1 fnum at bt ct dt et Fps Gps Hpf Mpf Gpf;

  % select the subproblem according to its utility.
  % if params.dynamic is not true, then no selection is used.
%   if (params.dynamic)
%     selindex = util_select();
%   else
    selindex = 1:length(subproblems);  
%   end
  
  %disp(selindex(1:5));

  global selectionSize;
  selectionSize = length(selindex);
  
  for i=1:length(selindex)
    index = selindex(i);
    
    %[r, rnduni]=crandom(rnduni);
    r = rand;
    updateneighbour = r < params.updateprob;
    %new point generation using genetic operations, and evaluate it.
    ind = genetic_op(index, updateneighbour, mop.domain);
    
    [obj,ind] = evaluate(mop, ind,itrCounter);
    %update the idealpoint.
    idealpoint = min(idealpoint, obj);
    %disp(evalCounter);
    %disp(obj);
    
    %updation.
    update(index, ind, updateneighbour);

    %clear!
    clear ind obj updateneighbour;
  end
end

% update the index's neighbour with the given individual.
% index is the subproblem's index in the main population.
% ind is the individual structure.
% updatenieghbour is a bool determine whether the neighbourhood of index, or the whole population should be updated.
% this procedure is also governed by a parameter from params: params.updatenb, which determine how many subproblem
% should be updated at most by this new individual.
function update(index, ind, updateneighbour)
  global subproblems idealpoint params window1 step1 fnum at bt ct dt et Fps Gps Hpf Mpf Gpf;
  
  % collect the updation index
  if (updateneighbour)
    updateindex = subproblems(index).neighbour;
  else
    updateindex = 1:length(subproblems);
  end
  
  updateindex = random_shuffle(updateindex);
  time=0;
  
  for i=1:length(updateindex)
      idx = updateindex(i);
      updateweight = subproblems(idx).weight;
      
      newobj=subobjective(updateweight, ind.objective,  idealpoint, 'te');
      old=subobjective(updateweight, subproblems(idx).curpoint.objective,  idealpoint, 'te');
      
      if (newobj<old)
         subproblems(idx).curpoint=ind;
         %disp(sprintf('UP -- ID: %d, Type: %d, T: %d, UI: %d -- %f', index-1, updateneighbour, time, idx-1, ind.parameter));
         %time = time+1; 
      end
      %if (time>=params.updatenb)
       %   return;
      %end
  end
  
  % do the comparision.
%   updateweight = [subproblems(updateindex).weight];
%   oops = [subproblems(updateindex).curpoint];
%   newobj=subobjective(updateweight, ind.objective,  idealpoint, 'te');
%   oldobj=subobjective(updateweight, [oops.objective], idealpoint, 'te' );
%   C = newobj < oldobj;
%   
%   % find the one need to be updated in betterIndex.
%   betterIndex = find(C);
%   updateNumber = params.updatenb;
%   if (length(betterIndex)>updateNumber)
%       randp = randperm(length(betterIndex));
%       randp = randp(1:updateNumber);
%       betterIndex = betterIndex(randp);
%   end
%   
%   % do the updation.
%   [subproblems(updateindex(betterIndex)).curpoint]= deal(ind);
%   clear C newobj oops oldobj;
end

function y =terminate()
    global params itrCounter;
    y = itrCounter>params.iteration;
end

function status()
    global evalCounter itrCounter selectionSize subproblems;
    %average utility
    averageutil = mean([subproblems.utility]);
    %if (~rem(itrCounter, 30))
    disp(sprintf('Itr:%d\tSel:%d\tEval:%d\tUtilM:%1.4f', ...
        itrCounter, selectionSize, evalCounter, averageutil));
    %end
end
