%matlab code for generating IGD value for any pareto-optimal front
%'front' is the generated front which is to be compared
%'name' is the problem name
%'dim' is the dimension of the problem
%no  is the number of points in the actual pareto front
%this is  for minimization

function [metric ] = igd3(F, name, dim)       %'front' is the set of non-dominated points in f-space(pareto optimal front)
N1=length(F(:,1));                                %metric is the igd value
m1=length(F(1,:));
                                                %F is a matrix,where
                                                %each row contains the
                                                %objective function values of a single particle along coloumn is objective number

         
    front=F;                                          
no=5000;
m = length(front(:,1));
N = length(front(1,:));

switch lower(name)
     case {'r2_dtlz2'}
         load R2_DTLZ2_M5.dat
         no=length(R2_DTLZ2_M5);
         PF=R2_DTLZ2_M5';
     case {'r2_dtlz3'}
         load R2_DTLZ3_M5.dat
         no=length(R2_DTLZ3_M5);
         PF=R2_DTLZ3_M5';
     case {'wfg1_m5'}
         load WFG1_M5.dat
         no=length(WFG1_M5);
         PF=WFG1_M5';
     otherwise 
         error('Undefined test problem name');
 end



%[PF, PS] = pareto(name, no, dim);
%no=length(PF);
value=0;
for i=1:no
    p=PF(:,i);
    tempdist=[];
    for j=1:N
        q = front(:,j);
        r=(p-q).^2;
        s=sum(r);
        s=sqrt(s);
        tempdist=[tempdist s];
    end 
  [Y M]=min(tempdist);
  value=value+Y;
end
value=value/no;
metric=value;
    
    