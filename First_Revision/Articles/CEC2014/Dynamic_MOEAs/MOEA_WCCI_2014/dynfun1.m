n=50;
t=1:n;
wt=t*2*pi/n;
g=sin(wt);
g=abs(g);
k=0.2+0.6*g;
x=0:0.01:1;
h=[1/5+24/5*g];
for i=1:n
f1=x+g(i);
f2=1-k(i).*(f1-g(i)).^h(i)+g(i);
if i<=15
    gg='r';
    plot(f1,f2,gg);
else
    gg='b';
    plot(f1,f2,gg);
end

hold on
end
hold off