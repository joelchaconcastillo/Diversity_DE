function mop = testmop( testname, dimension )
%Get test multi-objective problems from a given name. 
%   The method get testing or benchmark problems for Multi-Objective
%   Optimization. The test problem will be encapsulated in a structure,
%   which can be obtained by function get_structure('testmop'). 
%   User get the corresponding test problem, which is an instance of class
%   mop, by passing the problem name and optional dimension parameters.

    mop=get_structure('testmop');
    
    switch lower(testname)
        case {'fda1','zjz1','zdt11','fda2','fda3','zjz2','dyn1','dyn1_ps','dyn1_pf','ud1','ud2','ud3','ud4','ud5','ud6','ud7','ud8','ud9','ud10','ud11','ud12','ud13','ud14'}
            mop=dynmoprob(mop,testname,dimension);
        case {'cf1_var1','cf1_var2'}
            mop=dynconsprob(mop,testname,dimension);
        otherwise 
            error('Undefined test problem name');                
    end 
end


function p=dynmoprob(p,testname,dim)
p.name=testname;
p.pd=dim;
p.od=2;

switch testname
    case {'fda1','zdt11','fda2','zjz2','fda3'}
        bound1=[0 1];
        bound2=repmat([-1 1],dim-1,1);
        p.domain=[bound1;bound2];
%     case{'dyn1','dyn1_pf'}
%         bound1=[0 1];
%         bound2=repmat([-1 1],dim-1,1);
%         p.domain=[bound1;bound2];
%     case{'dyn1_ps'}
%         bound1=[0 1];
%         bound2=repmat([-2 2],dim-1,1);
%         p.domain=[bound1;bound2];
   case{'ud1','ud5','ud4','ud10','ud7','ud11'}
        bound1=[0 1];
        bound2=repmat([-1 1],dim-1,1);
        p.domain=[bound1;bound2];
  case{'ud6','ud12'}
        bound1=[0 1; 0 1];
        bound2=repmat([-3 3],dim-2,1);
        p.domain=[bound1;bound2];
        p.od=3;
  case{'ud2','ud8','ud13'}
        bound1=[0 1];
        bound2=repmat([-2 2],dim-1,1);
        p.domain=[bound1;bound2];
 case{'ud3','ud9','ud14'}
        bound1=[0 1];
        bound2=repmat([-1 2],dim-1,1);
        p.domain=[bound1;bound2];
    case {'zjz1'}
        bound1=[0 1];
        bound2=repmat([-1 2],dim-1,1);
        p.domain=[bound1;bound2];
end
p.func=benchmark_func(upper(testname));
end
function p=dynconsprob(p,testname,dim)
p.name=testname;
p.pd=dim;
p.od=2;
switch testname
    case {'cf1_var1','cf1_var2'}
        bound1=[0 1];
        bound2=repmat([-1 2],dim-1,1);
        p.domain=[bound1;bound2];
end
p.func=benchmark_func(upper(testname));
end
