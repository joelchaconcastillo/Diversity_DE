tao=(5*[0:20]);
% tao(length(tao))=[];
t=0:104;
j=1;
sine=[];

for i=1:max(t)+1
    if (mod(i-1,5)==0)
      ideal(i)=sin(tao(j)*pi/(10*5));
      sine=[sine ideal(i)];
      j=j+1;
    else 
      ideal(i)=ideal(i-1);
    end
end
idealtao=sin(tao*pi/(10));

% plot(t,ideal,'.-')
% hold on
% plot(tao,sine,'o-')
% hold off
opmean=0.05;

newpoint=[-0.3];
depoint=newpoint;
j=2;
for i=2:max(t)+1
    if (mod(i-1,5)==0)
        newpoint(j)=depoint(i-1)+(depoint(i-1)-newpoint(j-1));
        if newpoint(j)<= ideal(i)
            depoint(i)=newpoint(j)+opmean*abs(newpoint(j)-ideal(i));
        elseif newpoint(j)> ideal(i)
            depoint(i)=newpoint(j)-opmean*abs(newpoint(j)-ideal(i));
        end
        j=j+1;
    else
        if depoint(i-1)<= ideal(i)
            depoint(i)=depoint(i-1)+opmean*abs(depoint(i-1)-ideal(i));
        elseif depoint(i-1)> ideal(i)
            depoint(i)=depoint(i-1)-opmean*abs(depoint(i-1)-ideal(i));
        end
    end        
end

newpoint1=[-0.3];
depoint1=newpoint1;
j=2;
for i=2:length(t)
    
     if (mod(i-1,5)==0)&&j==2
        newpoint1(j)=depoint1(i-1)+(depoint1(i-1)-newpoint1(j-1));
        if newpoint(j)<= ideal(i)
            depoint1(i)=newpoint1(j)+opmean*abs(newpoint1(j)-ideal(i));
        elseif newpoint(j)> ideal(i)
            depoint1(i)=newpoint1(j)-opmean*abs(newpoint1(j)-ideal(i));
        end
        j=j+1;
        
     elseif (mod(i-1,5)==0)&&j>2
         
          x=abs(newpoint1(j-1)-newpoint1(j-2))-abs(depoint1(i-1)-newpoint1(j-1));
            if rand <0.5
                  epsilon=1+tanh(x);%+abs(randn);
            else
                  epsilon=1+(x)*abs(randn);
            end
            
            newpoint1(j)=depoint1(i-1)+epsilon*(depoint1(i-1)-newpoint1(j-1));   
        if newpoint1(j)<= ideal(i)
            depoint1(i)=newpoint1(j)+opmean*abs(newpoint1(j)-ideal(i));
        elseif newpoint1(j)> ideal(i)
            depoint1(i)=newpoint1(j)-opmean*abs(newpoint1(j)-ideal(i));
        end
        j=j+1;
        
     else
       
        if depoint1(i-1)<= ideal(i)
            depoint1(i)=depoint1(i-1)+opmean*abs(depoint1(i-1)-ideal(i));
        elseif depoint1(i-1)> ideal(i)
            depoint1(i)=depoint1(i-1)-opmean*abs(depoint1(i-1)-ideal(i));
        end
     end
end

plot(t,ideal,'.-',t,depoint,'*-',t,depoint1,'^-')
hold on
plot(tao,newpoint,'o',tao,newpoint1,'b+')
hold off

error=sum((ideal-depoint).^2)
error=sum((ideal-depoint1).^2)
% plot(t,ideal,'o-',t,depoint,'*-',t,depoint1,'^-');
