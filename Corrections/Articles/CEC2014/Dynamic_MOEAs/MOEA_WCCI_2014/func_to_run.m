%% variables to change at each run
%window in loadparams(change FEs and iteration),moead_basic,benchmark_func
%step in benchmark_func
%%
clc
global window1 step1 metric1 fnum at bt ct dt et Fps Gps Hpf Mpf Gpf;



dimension=10;
objectivedim=2;
name='ud14';  % use proper function name as given in benchmark_func.m where each function is defined
mop=testmop(name,dimension);
igdval=[];
metric=[];
num=300;
window1=5; step1=5; fnum=0; at=1; bt=1; ct=1; dt=1; et=1;  G1=0; G2=0; Hpf=0.5; Mpf=0.5; Gpf=0.0;  %use proper parameters for the functions
%initial setting for random ud13,ud14
Gps=G2;
F=ceil(dimension*G1);
Fps=F;
for i=1:1

[pareto metric1]=moead_pf(mop,{});% Use moead_ps, moead_pf, moead_newps, moead_newpf as 4 different algorithms
%[pareto ]=moead2(mop,{});
%front=feasible(pareto);
front=[pareto.objective];
%xvector=[pareto.parameter];
scatter(front(1,:),front(2,:),'o');

metric=[metric;  metric1];
avg=mean(metric1,2)
% sd=std(metric1,2);

end
