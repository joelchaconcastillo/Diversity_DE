%matlab code for generating IGD value for any pareto-optimal front
%'front' is the generated front which is to be compared
%'name' is the problem name
%'dim' is the dimension of the problem
%no  is the number of points in the actual pareto front
%this is  for minimization

function [metric ] = igd2(front, name, dim,iter)       %'front' is the set of non-dominated points in f-space(pareto optimal front)
global window1 step1  at bt ct dt et Fps Gps Hpf Mpf Gpf ;
% N1=length(F(:,1));                                %metric is the igd value
% m1=length(F(1,:));
                                                %F is a matrix,where
                                                %each row contains the
                                                %objective function values of a single particle along coloumn is objective number

         
                                         
no=5000;
% m = length(front(:,1));
N = length(front(1,:));

%switch lower(name)
 %    case {'r2_dtlz2'}
  %       load R2_DTLZ2_M5.dat
   %      no=length(R2_DTLZ2_M5);
    %     PF=R2_DTLZ2_M5';
    % case {'r2_dtlz3'}
     %    load R2_DTLZ3_M5.dat
      %   no=length(R2_DTLZ3_M5);
       %  PF=R2_DTLZ3_M5';
     %case {'wfg1_m5'}
      %   load WFG1_M5.dat
       %  no=length(WFG1_M5);
        % PF=WFG1_M5';
     %otherwise 
      %   error('Undefined test problem name');
 %end


switch name
    case {'fda1','zdt11','zjz1'}
         f1=linspace(0,1,no);
         f2=ones(1,no)-(f1.^0.5);
         PF=[f1;f2];
    case {'fda2','zjz2'}
        n=dim;
        f1=linspace(0,1,no);
        t=(floor(iter/window1))/step1;
        H=0.75+0.7*sin(0.5*3.14*t);
        Htemp=H*ones(n-1,1);
        k=(Htemp).^2;
        arbit=(sum(k)+H)^-1;
        
        f2=ones(1,no)-f1.^arbit;
        PF=[f1;f2];
    case{'fda3'}
        f1=linspace(0,1,no);
        t=(floor(iter/window1))/step1;
        G=abs(sin(0.5*3.14*t));
        f2=(1+G).*(ones(1,no)-(f1./(1+G)).^0.5);
        PF=[f1;f2];
   
    case {'cf1_var1','cf1_var2','dyn1','dyn1_pf','dyn1_ps','ud1','ud2','ud3','ud4','ud5','ud6','ud7','ud8','ud9','ud10','ud11','ud12','ud13','ud14'}
        PF=pareto(name,no,dim,iter);
end

    no=length(PF);

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
  clear M
  value=value+Y;
end
value=value/no;
metric=value;
    
    