%matlab code for generating Hypervolume metric value for any pareto-optimal front
%'front' is the generated front which is to be compared
%'name' is the problem name
%'od' is the dimension of the problem
%no  is the number of points in the actual pareto front
%this is  for minimization
function [indicator]=hvmetric(front,name,iter,od)
global window1 step1
h=1;
for i=1:od
    mx=max(front(i,:));
    mn=min(front(i,:));
    diff=mx-mn;
    h=h*diff;
end
href=1;
no=1000;
PF =pareto(name,no,od,iter);
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
end


for i=1:od
    mx=max(PF(i,:));
    mn=min(PF(i,:));
    diff=mx-mn;
    href=href*diff;
end
indicator=href-h;
end