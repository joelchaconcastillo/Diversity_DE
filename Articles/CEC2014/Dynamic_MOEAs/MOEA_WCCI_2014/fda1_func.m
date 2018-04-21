function [f c]=fda1_func(x,iter)
x=x';
window=5; %time window for which the process remains constant
step=5;  %controls the distance between 2 paretos nT
n=length(x);
f1=x(1);
t=(floor(iter/window))/step;
G=sin(0.5*3.14*t);
temp=x(2:n);
Gtemp=G*ones(n-1,1);
k=(temp-Gtemp).^2;
arbit=sum(k);
g=1+arbit;
f2=1-(f1/g)^0.5;
f=[f1+f2];
c=[0 0];
end