
f1=0:0.01:1;
n=10;
t=1:n;
wt=t*2*pi/n;
g=sin(wt);
h=0.5*(5/2+3/2*g);
for i=1:n
f2mat(i,:)=1-f1.^h(i);
plot(f1,f2mat(i,:));
hold on
end
hold off