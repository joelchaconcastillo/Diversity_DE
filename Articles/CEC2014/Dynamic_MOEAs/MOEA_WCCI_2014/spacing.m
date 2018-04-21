function [metric] = spacing(front)
N=length(front(1,:));           %M is the no of objectives and N popsize 
M=length(front(:,1));
dist=[];
for i=1:N
    tempdist=[];
    for j=1:N
        tt=[];
      for k=1:M
        t=abs(front(k,i)-front(k,j));
        
        tt=[tt;t ];
      end
        temp= sum(tt);
        if i==j
            temp=Inf;
        end
      tempdist=[tempdist temp];
    end
      val=min(tempdist);    
     dist=[dist;val];
end
dist;
d=mean(dist);
dist=sum((dist-repmat(d,N,1)).^2)/N;
metric=sqrt(dist);
end