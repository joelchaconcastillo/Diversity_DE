load nps5
a=[mean(metric(1,:)) mean(metric(2,:))]
x=mean(a)
y=std(a)
fprintf('%d %d',x,y)