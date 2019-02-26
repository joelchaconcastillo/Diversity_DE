%==========================================================================
%    MVMO-based solution for CEC 2016 Special Session and Competition  
%    on Learning-based Real-Parameter Single Objective Optimization
%==========================================================================

%                                  Reference
%--------------------------------------------------------------------------
% [1] J. J. Liang, B. Y. Qu, P. N. Suganthan, 
%     "Problem Definitions and Evaluation Criteri for the CEC 2015 
%      Competition on Learning-based Real-Parameter Single Objective 
%      Optimization," Nov. 2014
% [Online] Available at http://www.ntu.edu.sg/home/epnsugan/
%--------------------------------------------------------------------------

%By: Dr.-Ing. José L. Rueda & Prof. István Erlich
%Email: j.l.ruedatorres@tudelft.nl, istvan.erlich@uni-due.de
%Date: 12.06.2016 


delete(gcp)


close all
clear all
clc


global proc


proc.fbest_stats=zeros(15,5);
proc.mvmo_time=zeros(1,1);
proc.mvmo_time1=zeros(1,1);

%For algorithm Complexity
tStart1 = cputime;
for i=1:1000000
    x= 0.55d0+i;
    x=x+x; 
    x=x/2;
    x=x*x; 
    x=sqrt(x); 
    x=log(x);
    x=exp(x); 
    y=x/(x+2);
end
T0 = cputime-tStart1;

ranges = [-100,100];

%                Preliminary definitions & Parallelization
%--------------------------------------------------------------------------
refresh=1000; %Printing step
algorithm_name='mvmo';
algorithm_hd=str2func(algorithm_name);
test_bed_OPF_hd=str2func('test_func_calc');
args{1}=refresh;
args{2}=algorithm_name;
args{3}=10; %30; %50; %100; %Problem dimension
args{4}=10000*args{3}; %Maximum function evaluations
% args{5}=repmat(-100,1,args{3}); %Min limits of control variables 
% args{6}=repmat(100,1,args{3});%Max limits of control variables  
args{7}=[100:100:1500]'; %Theoretical optimum value
args{8}=51; %Optimization runs
args{9}=1;%Printing 1=yes, 0=no
 
run_in_parallel=1; %1 = yes,  0 = no

NumWorkers=4;


if run_in_parallel
    c = parcluster;
    c.NumWorkers=NumWorkers;
    if run_in_parallel
        parpool(c);
    end
end
%--------------------------------------------------------------------------
%                           Running optimization
%--------------------------------------------------------------------------
fnum_count=zeros(15,1);
icont=1;


for func_num=1:15 %From 1 to 15 - Problem number
    
     
    lu_lim = ranges;
    args{5}=lu_lim(1)*ones(1,args{3}); %Min limits of control variables
    args{6}=lu_lim(2)*ones(1,args{3}); %Max limits of control variables 
        
    test_bed_OPF_hd(func_num,1,icont,args); 
    %For algorithm Complexity
    tStart2 = cputime;
    parfor krun=1:args{8}
        test_bed_OPF_hd(func_num,0,icont,args);
        % Call to your implementation.
        feval(algorithm_hd,test_bed_OPF_hd,func_num,krun,icont,args);
%         fprintf('Run %d finished.\n',krun);
    end
    proc.mvmo_time = (cputime-tStart2)/args{8};
    test_bed_OPF_hd(func_num,987,icont,args);   
       
    
    fnum_count(func_num,:)=func_num;  
    
end


%--------------------------------------------------------------------------


% Parallelization
% Deactivation of 
% shared-memory session.
if run_in_parallel
    poolobj = gcp('nocreate');
    delete(poolobj);
end


tot_score=sum(proc.fbest_stats(:,3))+sum(proc.fbest_stats(:,4));
fbest_stats=[fnum_count, proc.fbest_stats];
mvmo_time=[T0,proc.mvmo_time1,proc.mvmo_time,(proc.mvmo_time-proc.mvmo_time1)/T0];
save(proc.problem_filename_stats, 'fbest_stats','-ASCII')
save(proc.problem_filename_time, 'mvmo_time','-ASCII')

