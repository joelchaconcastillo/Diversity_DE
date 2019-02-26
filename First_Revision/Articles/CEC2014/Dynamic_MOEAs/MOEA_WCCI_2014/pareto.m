% pareto.m
% 
% The Matlab source codes to generate the PF and the PS of the test
%   instances for CEC 2009 Multiobjective Optimization Competition.
% 
% Usage: [pf, ps] = pareto(problem_name, no_of_points, variable_dim)
% 
% Please refer to the report for more information.
% 
% History:
%   v1 Sept.08 2008

function [pf, ps] = pareto(name, no, dim,iter)
global window1 step1 at bt ct dt et Fps Gps Hpf Mpf Gpf;
    if nargin<3, dim = 3; end
    if nargin<2, no  = 500; end
    switch name
         
        case {'cf1_var1'}
            window1=5;
            step1=5;
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
            Ftemp=1.5+G;
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = (0:1:20)/20.0;
            pf(2,:)     = 1-pf(1,:).^Ftemp;
            ps          = [];  
        case {'cf1_var2'}
            window1=5;
            step1=5;
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
            
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = (0:1:20)/20.0;
            pf(2,:)     = 1+abs(G)-pf(1,:);
            ps          = [];  
        case {'dyn1'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
            N=5;
            H=1/N+(N-1/N)*abs(G);
            K=0.2+0.6*abs(G);
            pf          = zeros(2,no);
            tempo=linspace(0,1,no);
            pf(1,:)     = tempo+abs(G);
            pf(2,:)     = 1+abs(G)-(tempo).^H;
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
         case {'dyn1_ps'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
            N=2;
%             H=1/N+(N-1/N)*abs(G);
            H=0.5+abs(G);
            K=0.2+0.6*abs(G);
            pf          = zeros(2,no);
            tempo=linspace(0,1,no);
            pf(1,:)     = tempo;
            pf(2,:)     = 1-(tempo).^0.5;
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
          case {'ud1','ud2','ud3'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
%             N=2;
%             H=1/N+(N-1/N)*abs(G);
%             H=0.5+abs(G);
%             K=0.5+abs(G);
            pf          = zeros(2,no);
            tempo=linspace(0,1,no);
            pf(1,:)     = tempo+abs(G);
            pf(2,:)     = 1+abs(G).^1-(tempo);
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
       case {'ud4'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
           
%             H=1/N+(N-1/N)*abs(G);
%             H=0.5+abs(G);
%             K=0.5+abs(G);
%             pf          = zeros(2,no);
%             tempo=linspace(0,1,no);
% %             pf(1,:)     = tempo+abs(G);
% %             pf(2,:)     = 1+abs(G)-K*(tempo).^H;
%             pf(1,:)     = tempo;
%             pf(2,:)     = 1-tempo;
%             ps          = zeros(dim,no);
%             ps(1,:)     = linspace(0,1,no);
%             ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = ((0:1:20)/20.0)+G;
            pf(2,:)     = 1-pf(1,:)+2*G;
            ps          = zeros(dim,no);
            ps(1,:)     = pf(1,:);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]));
        case 'ud5'
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
            G=abs(G);
           
%             H=1/N+(N-1/N)*abs(G);
            H=0.5+abs(G);
            K=0.5+abs(G);
            num                     = floor(no/3);
            pf                      = zeros(2,no);
            pf(1,1:num)             = 0.0;
            pf(1,(num+1):(2*num))   = linspace(0.25,0.5,num);
            pf(1,(2*num+1):no)      = linspace(0.75,1.0,no-2*num);
            pf(1,:)                 = pf(1,:)+G; 
            pf(2,:)                 = 1-pf(1,:)+2*G;
            
            ps                      = zeros(dim,no);
            ps(1,:)                 = pf(1,:);
            ps(2:dim,:)             = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
        case {'ud6'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*pi*t);
            F=1+abs(G);
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = F.*(cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:)));
            pf(2,:)     = F.*(cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:)));
            pf(3,:)     = F.*(sin(0.5*pi*ps(1,:)));   
            clear s t;
          case {'ud7','ud8','ud9'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
           % N=2;
%             H=1/N+(N-1/N)*abs(G);
            %H=0.5+abs(G);
            K=0.5+abs(G);
            pf          = zeros(2,no);
            tempo=linspace(0,1,no);
            pf(1,:)     = tempo;
            pf(2,:)     = 1-K.*(tempo).^K;
            ps          = zeros(dim,no);
            ps(1,:)     = linspace(0,1,no);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
       case {'ud10'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*pi*t);
           
%             H=1/N+(N-1/N)*abs(G);
            H=0.5+abs(G);
%             K=0.5+abs(G);
%             pf          = zeros(2,no);
%             tempo=linspace(0,1,no);
% %             pf(1,:)     = tempo+abs(G);
% %             pf(2,:)     = 1+abs(G)-K*(tempo).^H;
%             pf(1,:)     = tempo;
%             pf(2,:)     = 1-tempo;
%             ps          = zeros(dim,no);
%             ps(1,:)     = linspace(0,1,no);
%             ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
            no          = 21;
            pf          = zeros(2,no);
            pf(1,:)     = ((0:1:20)/20.0)+abs(G);
            pf(2,:)     = 1-pf(1,:).*H+2*(G);
            ps          = zeros(dim,no);
            ps(1,:)     = pf(1,:);
            ps(2:dim,:) = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]));
    
    
    case 'ud11'
            t=(floor(iter/window1))/step1;
            G=sin(0.5*3.14*t);
            G=abs(G);
           
%             H=1/N+(N-1/N)*abs(G);
            %H=0.5+abs(G);
            K=0.5+abs(G);
            num                     = floor(no/3);
            pf                      = zeros(2,no);
            pf(1,1:num)             = 0.0;
            pf(1,(num+1):(2*num))   = linspace(0.25,0.5,num);
            pf(1,(2*num+1):no)      = linspace(0.75,1.0,no-2*num);
            pf(1,:)                 = pf(1,:)+G; 
            pf(2,:)                 = 1-pf(1,:).*K+2*G;
            
            ps                      = zeros(dim,no);
            ps(1,:)                 = pf(1,:);
            ps(2:dim,:)             = sin(6.0*pi*repmat(ps(1,:),[dim-1,1]) + repmat((2:dim)',[1,no])*pi/dim);
   
    
     case {'ud12'}
            t=(floor(iter/window1))/step1;
            G=sin(0.5*pi*t);
            F=1+abs(G);
            num         = floor(sqrt(no));
            no          = num*num;
            [s,t]       = meshgrid(linspace(0,1,num),linspace(0,1,num));
            ps          = zeros(dim,no);
            ps(1,:)     = reshape(s,[1,no]);
            ps(2,:)     = reshape(t,[1,no]);            
            ps(3:dim,:) = 2.0*repmat(ps(2,:),[dim-2,1]).*sin(2.0*pi*repmat(ps(1,:),[dim-2,1]) + repmat((3:dim)',[1,no])*pi/dim);             
            pf          = zeros(3,no);
            pf(1,:)     = F.*(cos(0.5*pi*ps(1,:)).*cos(0.5*pi*ps(2,:)))+abs(G);
            pf(2,:)     = F.*(cos(0.5*pi*ps(1,:)).*sin(0.5*pi*ps(2,:)))+abs(G);
            pf(3,:)     = F.*(sin(0.5*pi*ps(1,:)))+abs(G);   
            clear s t;

    case {'ud13','ud14'}
            %t=(floor(iter/window1))/step1;
            %G=sin(0.5*3.14*t);
%             N=2;
%             H=1/N+(N-1/N)*abs(G);
%             H=0.5+abs(G);
%             K=0.5+abs(G);

            pf          = zeros(2,no);
            tempo=linspace(0,1,no);
            pf(1,:)     = tempo+Gpf;
            pf(2,:)     = 1+Gpf-Mpf.*(tempo).^Hpf;
    end
    
    
end
