function fobj= benchmark_func(testname) %x is a coloumn vector
%global window1 step1;
switch lower(testname)
    case {'fda1'}
        fobj=@fda1_func;
     case {'fda2'}
        fobj=@fda2_func; 
     case {'fda3'}
        fobj=@fda3_func; 
    case {'cf1_var1'}
        fobj=@cf1_var1_func;
    case {'cf1_var2'}
        fobj=@cf1_var2_func;
    case {'ud1'}
        fobj=@ud1func;
    case {'ud2'}
        fobj=@ud2func;
    case {'ud3'}
        fobj=@ud3func;
    case {'ud4'}
        fobj=@ud4func;
    case {'ud5'}
        fobj=@ud5func;
   case {'ud6'}
        fobj=@ud6func;
    case {'ud7'}
        fobj=@ud7func;
   case {'ud8'}
        fobj=@ud8func;
   case {'ud9'}
        fobj=@ud9func; 
   case {'ud10'}
        fobj=@ud10func;   
   case {'ud11'}
        fobj=@ud11func;  
   case {'ud12'}
        fobj=@ud12func;      
    case {'ud13'}
        fobj=@ud13func;
    case {'ud14'}
        fobj=@ud14func;    
        
    otherwise 
        error('Undefined test problem name');                
end 





end


function f=fda1_func(x,iter)
global window1 step1;

% window1=5;
% step1=5;
% window1=5; %time window1 for which the process remains constant
% step1=5;  %controls the distance between 2 paretos nT
n=length(x);
f1=x(1);
t=(floor(iter/window1))/step1;
G=sin(0.5*3.14*t);
temp=x(2:n);
Gtemp=G*ones(n-1,1);
k=(temp-Gtemp).^2;
arbit=sum(k);
g=1+arbit;
f2=g*(1-(f1/g)^0.5);
f=[f1;f2];
end

function f=fda2_func(x,iter)
global window1 step1;

% window1=5; %time window1 for which the process remains constant
% step1=5;  %controls the distance between 2 paretos nT
n=length(x);
f1=x(1);
t=(floor(iter/window1))/step1;
H=0.75+0.7*sin(0.5*3.14*t);
temp=x(2:n);
g=1+sum(temp.^2);
Htemp=H*ones(n-1,1);
k=(temp-Htemp).^2;
arbit=(sum(k)+H)^-1;
f2=g*(1-(f1/g)^arbit);
f=[f1;f2];
end

function f=fda3_func(x,iter)
global window1 step1;
% window1=5; %time window1 for which the process remains constant
% step1=5;  %controls the distance between 2 paretos nT
n=length(x);
t=(floor(iter/window1))/step1;
G=abs(sin(0.5*3.14*t));
F=2*G;%10^(2*sin(0.5*3.14*t))
f1=x(1)^F;

Gtemp=G*ones(n-1,1);
temp=x(2:n);
g=1+G+sum((temp-Gtemp).^2);
f2=g*(1-(f1/g)^0.5);
f=[f1;f2];
end
function [f,c]=cf1_var1_func(x,iter)
global window1 step1;
    
    [dim, num]   = size(x);
     window1=5; %time window1 for which the process remains constant
    step1=5;  %controls the distance between 2 paretos nT
    t=(floor(iter/window1))/step1;
    G=sin(0.5*3.14*t);
    Gtemp=G*ones(dim-1,num);
    Ftemp=1.5*ones(1,num)+G*ones(1,num);
    a            = 1.0;
    N            = 10.0;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:).^Ftemp + 2.0*tmp2/size(2:2:dim,2);
    c(1,:)       = y(1,:).^Ftemp + y(2,:) - a*abs(sin(N*pi*(y(1,:).^Ftemp-y(2,:)+1.0))) - 1.0;
    f=[y(1,:);y(2,:)];
end
function [f,c]=cf1_var2_func(x,iter)
global window1 step1;
    
    [dim, num]   = size(x);
     window1=5; %time window1 for which the process remains constant
    step1=5;  %controls the distance between 2 paretos nT
    t=(floor(iter/window1))/step1;
    G=sin(0.5*3.14*t);
    Gtemp=G*ones(dim-1,num);
    %Ftemp=1.5*ones(1,num)+G*ones(1,num);
    a            = 1.0;
    N            = 10.0;
    [dim, num]   = size(x);
    Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(0.5+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)       = x(1,:)       + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)       = 1.0+abs(Gtemp(1,:)) - x(1,:) + 2.0*tmp2/size(2:2:dim,2);
    c(1,:)       = y(1,:) + y(2,:) - a*abs(sin(N*pi*(y(1,:).^Ftemp-y(2,:)+1.0))) - 1.0-abs(Gtemp(1,:));
    f=[y(1,:);y(2,:)];
end



function f = ud1func(x,iter)
global window1 step1;
    [dim, num]  = size(x);

    t=(floor(iter/window1))/step1;

    G=sin(0.5*pi*t);
    Gtemp=G*ones(1,num);
    F=ceil(dim*G);
    Ftemp=F*ones(dim-1,num);

    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)+abs(Gtemp) + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - (x(1,:)).^1+abs(Gtemp) + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
end
function f = ud2func(x,iter)
global window1 step1;
    [dim, num]  = size(x);
    t=(floor(iter/window1))/step1;
    G=sin(0.5*pi*t);
    Gtemp=G*ones(1,num);
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)+abs(Gtemp) + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - (x(1,:).^1)+abs(Gtemp) + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
end
function f = ud3func(x,iter)
global window1 step1;
    [dim, num]  = size(x);

    t=(floor(iter/window1))/step1;

    G=sin(0.5*pi*t);
    Gtemp=G*ones(1,num);


    Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(1+Gtemp+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)      = x(1,:)+abs(Gtemp) + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - (x(1,:)).^1+abs(Gtemp) + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
end

function f = ud4func(x,iter)
global window1 step1;
[dim, num]  = size(x);
     t=(floor(iter/window1))/step1;
     G=sin(0.5*pi*t);
     Gtemp=G*ones(1,num);
    N           = 10.0;
    E           = 0.1;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(2:dim,:)  = x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    H           = zeros(dim,num);
    H(2:dim,:)  = 2.0*Y(2:dim,:).^2 - cos(4.0*pi*Y(2:dim,:)) + 1.0;
    tmp1        = sum(H(3:2:dim,:));  % odd index
    tmp2        = sum(H(2:2:dim,:));  % even index
    tmp         = (0.5/N+E)*abs(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp));
    y(1,:)      = x(1,:)      + tmp + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - x(1,:)+ tmp + 2.0*tmp2/size(2:2:dim,2);
    f=[y(1,:);y(2,:)];
end
function y = ud5func(x,iter)
global window1 step1;
[dim, num]   = size(x);
t=(floor(iter/window1))/step1;
     G=sin(0.5*pi*t);
     Gtemp=G*ones(1,num);
    N            = 2.0;
    E            = 0.1;
    Y            = zeros(dim,num);
    Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1         = zeros(dim,num);
    tmp1(2:dim,:)= Y(2:dim,:).^2;
    tmp2         = zeros(dim,num);
    tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
    tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
    tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
    tmp          = max(0,(1.0/N+2.0*E)*(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp)));
    y(1,:)       = x(1,:)       + tmp + 2.0*tmp11/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:) + tmp + 2.0*tmp21/size(2:2:dim,2);
    clear Y tmp1 tmp2;
    f=[y(1,:);y(2,:)];
end
function y = ud6func(x,iter)
    global window1 step1;
[dim, num]   = size(x);
t=(floor(iter/window1))/step1;
     G=sin(0.5*pi*t);
     F=1+abs(G);
     Gtemp=G*ones(1,num);
     Ftemp=F*ones(1,num);
   
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:)-Gtemp - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) + 2.0*tmp1/size(4:3:dim,2));
    y(2,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) +2.0*tmp2/size(5:3:dim,2));
    y(3,:)      = Ftemp.*(sin(0.5*pi*x(1,:))                     + 2.0*tmp3/size(3:3:dim,2));
    clear Y;
end
function f = ud7func(x,iter)
global window1 step1;
    [dim, num]  = size(x);
%     tmp         = zeros(dim,num);
    t=(floor(iter/window1))/step1;
%     N=2;
    G=sin(0.5*pi*t);
%    Gtemp=G*ones(1,num);
    F=ceil(dim*G);
    Ftemp=F*ones(dim-1,num);
%     K=0.5+abs(G);
%     Ktemp=K*ones(1,num);
%     H=1/N+(N-1/N)*abs(G);
     H=0.5+abs(G);
%     Htemp=H*ones(1,num);
    
%     tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(F+repmat((2:dim)',[1,num])))-G).^2;
%     tmp1        = sum(tmp(3:2:dim,:));  % odd index
%     tmp2        = sum(tmp(2:2:dim,:));  % even index
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)+ + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - H.*(x(1,:)).^H + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
end
function f = ud8func(x,iter)
global window1 step1;
    [dim, num]  = size(x);
    t=(floor(iter/window1))/step1;
    G=sin(0.5*pi*t);
    Gtemp=G*ones(1,num);
    %F=ceil(dim*G);
    %Ftemp=F*ones(dim-1,num);
    K=0.5+abs(G);
    %Ktemp=K*ones(1,num);
%     H=1/N+(N-1/N)*abs(G);
    %H=0.5+abs(G);
    %Htemp=H*ones(1,num);
    
%     tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(F+repmat((2:dim)',[1,num])))-G).^2;
%     tmp1        = sum(tmp(3:2:dim,:));  % odd index
%     tmp2        = sum(tmp(2:2:dim,:));  % even index
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(repmat((2:dim)',[1,num])))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)+ 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - K.*(x(1,:)).^K + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
end
function f = ud9func(x,iter)
global window1 step1;
    [dim, num]  = size(x);
%     tmp         = zeros(dim,num);
    t=(floor(iter/window1))/step1;
%     N=2;
    G=sin(0.5*pi*t);
    Gtemp=G*ones(1,num);
%     F=ceil(dim*G);
%     Ftemp=F*ones(dim-1,num);
    K=0.5+abs(G);
%     Ktemp=K*ones(1,num);
%     H=1/N+(N-1/N)*abs(G);
%     H=0.5+abs(G);
%     Htemp=H*ones(1,num);
    
%     tmp(2:dim,:)= (x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(F+repmat((2:dim)',[1,num])))-G).^2;
%     tmp1        = sum(tmp(3:2:dim,:));  % odd index
%     tmp2        = sum(tmp(2:2:dim,:));  % even index
%     tmp         = zeros(dim,num);
%     tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
%     tmp1        = sum(tmp(3:2:dim,:));  % odd index
%     tmp(2:dim,:)= (x(2:dim,:)-Gtemp - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Ftemp+repmat((2:dim)',[1,num])))).^2;
%     tmp2        = sum(tmp(2:2:dim,:));  % even index
    Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(1+Gtemp+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gtemp).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)      = x(1,:) + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - K.*(x(1,:)).^K+ 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
end
function f = ud10func(x,iter)
global window1 step1;
[dim, num]  = size(x);
     t=(floor(iter/window1))/step1;
     G=sin(0.5*pi*t);
     K=0.5+abs(G);
     Gtemp=G*ones(1,num);
    N           = 10.0;
    E           = 0.1;
    [dim, num]  = size(x);
    Y           = zeros(dim,num);
    Y(2:dim,:)  = x(2:dim,:)- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    H           = zeros(dim,num);
    H(2:dim,:)  = 2.0*Y(2:dim,:).^2 - cos(4.0*pi*Y(2:dim,:)) + 1.0;
    tmp1        = sum(H(3:2:dim,:));  % odd index
    tmp2        = sum(H(2:2:dim,:));  % even index
    tmp         = (0.5/N+E)*abs(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp));
    y(1,:)      = x(1,:)      + tmp + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1- x(1,:).*K+ tmp + 2.0*tmp2/size(2:2:dim,2);
    f=[y(1,:);y(2,:)];
end
function y = ud11func(x,iter)
global window1 step1;
[dim, num]   = size(x);
t=(floor(iter/window1))/step1;
     G=sin(0.5*pi*t);
     K=0.5+abs(G);
     Gtemp=G*ones(1,num);
    N            = 2.0;
    E            = 0.1;
    Y            = zeros(dim,num);
    Y(2:dim,:)  = x(2:dim,:) - sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*repmat((2:dim)',[1,num]));
    tmp1         = zeros(dim,num);
    tmp1(2:dim,:)= Y(2:dim,:).^2;
    tmp2         = zeros(dim,num);
    tmp2(2:dim,:)= cos(20.0*pi*Y(2:dim,:)./sqrt(repmat((2:dim)',[1,num])));
    tmp11        = 4.0*sum(tmp1(3:2:dim,:)) - 2.0*prod(tmp2(3:2:dim,:)) + 2.0;  % odd index
    tmp21        = 4.0*sum(tmp1(2:2:dim,:)) - 2.0*prod(tmp2(2:2:dim,:)) + 2.0;  % even index
    tmp          = max(0,(1.0/N+2.0*E)*(sin(2.0*N*pi*x(1,:))-abs(2*N*Gtemp)));
    y(1,:)       = x(1,:)       + tmp + 2.0*tmp11/size(3:2:dim,2);
    y(2,:)       = 1.0 - x(1,:).*K + tmp + 2.0*tmp21/size(2:2:dim,2);
    clear Y tmp1 tmp2;
 %  f=[y(1,:);y(2,:)];
end
function y = ud12func(x,iter)
    global window1 step1;
[dim, num]   = size(x);
t=(floor(iter/window1))/step1;
     G=sin(0.5*pi*t);
     F=1+abs(G);
     Gtemp=G*ones(1,num);
     Ftemp=F*ones(1,num);
   
    Y           = zeros(dim,num);
    Y(3:dim,:)  = (x(3:dim,:)-Gtemp - 2.0*repmat(x(2,:),[dim-2,1]).*sin(2.0*pi*repmat(x(1,:),[dim-2,1]) + pi/dim*repmat((3:dim)',[1,num]))).^2;
    tmp1        = sum(Y(4:3:dim,:));  % j-1 = 3*k
    tmp2        = sum(Y(5:3:dim,:));  % j-2 = 3*k
    tmp3        = sum(Y(3:3:dim,:));  % j-0 = 3*k
    y(1,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*cos(0.5*pi*x(2,:)) +abs(Gtemp)+ 2.0*tmp1/size(4:3:dim,2));
    y(2,:)      = Ftemp.*(cos(0.5*pi*x(1,:)).*sin(0.5*pi*x(2,:)) + abs(Gtemp)+2.0*tmp2/size(5:3:dim,2));
    y(3,:)      = Ftemp.*(sin(0.5*pi*x(1,:))                     +abs(Gtemp)+ 2.0*tmp3/size(3:3:dim,2));
    clear Y;
end

function f = ud13func(x,iter)
global window1 step1 fnum at bt ct dt et Fps Gps Hpf Mpf Gpf;
    [dim, num]  = size(x);

    
    
if mod(iter,window1)==0 
Hpf;
end
    tmp         = zeros(dim,num);
    tmp(2:dim,:)= (x(2:dim,:) -Gps- sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Fps+repmat((2:dim)',[1,num])))).^2;
    tmp1        = sum(tmp(3:2:dim,:));  % odd index
    tmp(2:dim,:)= (x(2:dim,:)- Gps-sin(6.0*pi*repmat(x(1,:),[dim-1,1]) + pi/dim*(Fps+repmat((2:dim)',[1,num])))).^2;
    tmp2        = sum(tmp(2:2:dim,:));  % even index
    y(1,:)      = x(1,:)+abs(Gpf) + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - Mpf*(x(1,:)).^Hpf+abs(Gpf) + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];

end
function f = ud14func(x,iter)
global window1 step1 fnum at bt ct dt et Fps Gps Hpf Mpf Gpf;
    [dim, num]  = size(x);

    
    
if mod(iter,window1)==0
Hpf;
end
Y            = zeros(dim,num);
    Y(2:dim,:)   = (x(2:dim,:) - repmat(x(1,:),[dim-1,1]).^(1+Gps+1.5*(repmat((2:dim)',[1,num])-2.0)/(dim-2.0))-Gps).^2;
    tmp1         = sum(Y(3:2:dim,:));% odd index
    tmp2         = sum(Y(2:2:dim,:));% even index 
    y(1,:)      = x(1,:)+abs(Gpf) + 2.0*tmp1/size(3:2:dim,2);
    y(2,:)      = 1.0 - Mpf*(x(1,:)).^Hpf+abs(Gpf) + 2.0*tmp2/size(2:2:dim,2);
    clear tmp;
    f=[y(1,:);y(2,:)];
    
end