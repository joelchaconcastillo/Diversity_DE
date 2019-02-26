function metric=exmetric(front,set, name, dim,iter) 

global window1 step1;

no=1000;
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
   
    case {'cf1_var1','cf1_var2','dyn1','dyn1_pf','dyn1_ps','ud1','ud2','ud3','ud4','ud5','ud6','ud7','ud8','ud9','ud10','ud11','ud12'}
        PF=pareto(name,no,dim,iter);
end
    no=length(PF);
    [od N] = size(front);
    U=[];
   ind=[];
    for i=1:od
        [uu ff]=min(front(i,:));
        U=[U uu];
        ind=[ind ff];
    end
     Mv=[];
    for i=1:od
        t=set(:,ind(i));
       Mv=[Mv t'];
    end
    Mmat=[];
    for i=1:od
        for j=1:od
            if i==j
                Mmat(i,j)=U(i);
            else
                Mmat(i,j)=front(j,ind(i));
            end
        end
    end
    nadir=[];
    for i=1:od
        bb=max(Mmat(i,:));
        nadir=[nadir bb];
    end    
    val=norm(U-nadir);
    metric=0;
    for j=1:N
        cc=[];
        for i=1:no
            cc=[cc norm((PF(:,i)-front(:,j))/val)];
        end
        temp=min(cc);
        metric=metric+temp;
    end
    metric=metric/N;
end
    