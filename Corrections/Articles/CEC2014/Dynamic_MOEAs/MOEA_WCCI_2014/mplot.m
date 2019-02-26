function mplot(ps,name)
%X=[ps.objective];
X=ps;
[y i]=sort(X(1,:),'ascend');
z=X(2,:);
z=z(i);
%plot(y,z,'o')
%hold on
X1=pareto(name,5000,2);
[y1 i]=sort(X1(1,:),'ascend');
z1=X1(2,:);
z1=z1(i);
%plot(y,z,'-')
%hold off
plot(y,z,'go',y1,z1,'r-')
end