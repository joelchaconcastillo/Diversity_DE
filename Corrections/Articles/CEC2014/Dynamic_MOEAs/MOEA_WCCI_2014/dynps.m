n=10;
t=1:n;
wt=t*2*pi/n;
g=sin(wt);
f=ceil(30*g);
x1=0:.01:1;
for i=1:n
x2=sin(6*pi*x1+(pi/30)*(2+f(i)))-g(i);
if i<=5
plot(x1,x2,'c')
else
plot(x1,x2,'b')    
end
hold on
end
hold off