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

depoint=[-0.3];
storepoint=depoint;
j=1;
for i=2:max(t)+1
    if (mod(i-1,5)==0)&&j>2
        j=j+1;
        newpoint=storepoint(j-1)+(storepoint(j-1)-storepoint(j-2))
        if newpoint<= ideal(i)
            depoint(i)=newpoint+opmean*abs(newpoint-ideal(i));
        elseif newpoint> ideal(i)
            depoint(i)=newpoint-opmean*abs(newpoint-ideal(i));
        end
        
    else
        if depoint(i-1)<= ideal(i)
            depoint(i)=depoint(i-1)+opmean*abs(depoint(i-1)-ideal(i));
        elseif depoint(i-1)> ideal(i)
            depoint(i)=depoint(i-1)-opmean*abs(depoint(i-1)-ideal(i));
        end
        if (mod(i,5)==0)
            j=j+1;
            storepoint(j)=depoint(i);
        end
    end        
end

% % depoint1=[-0.3];
% % storepoint1=depoint1;
% % j=1;
% % for i=2:length(t)
% %     
% %      if (mod(i-1,5)==0)&&j==3
% %         j=j+1;
% %         newpoint1=storepoint1(j-1)+(storepoint1(j-1)-storepoint1(j-2));
% %         if newpoint1<= ideal(i)
% %             depoint1(i)=newpoint1+opmean*abs(newpoint1-ideal(i));
% %         elseif newpoint1> ideal(i)
% %             depoint1(i)=newpoint1-opmean*abs(newpoint1-ideal(i));
% %         end
% %         
% %         
% %      elseif (mod(i-1,5)==0)&&j>3
% %          j=j+1;
% %           x=abs(storepoint1(j-1)-storepoint1(j-2))-abs(storepoint1(j-3)-storepoint1(j-3));
% %             if rand <0.5
% %                   epsilon=1+tanh(x);%+abs(randn);
% %             else
% %                   epsilon=1+(x)*abs(randn);
% %             end
% %             
% %             newpoint1=storepoint1(j-1)+epsilon*(storepoint1(j-1)-storepoint1(j-2));   
% %         if newpoint1<= ideal(i)
% %             depoint1(i)=newpoint1+opmean*abs(newpoint1-ideal(i));
% %         elseif newpoint1> ideal(i)
% %             depoint1(i)=newpoint1-opmean*abs(newpoint1-ideal(i));
% %         end
% %         
% %         
% %      else
% %        
% %         if depoint1(i-1)<= ideal(i)
% %             depoint1(i)=depoint1(i-1)+opmean*abs(depoint1(i-1)-ideal(i));
% %         elseif depoint1(i-1)> ideal(i)
% %             depoint1(i)=depoint1(i-1)-opmean*abs(depoint1(i-1)-ideal(i));
% %         end
% %         if (mod(i,5)==0)
% %             j=j+1;
% %             storepoint1(j)=depoint1(i);
% %         end
% %      end
% % end

plot(t,ideal,'.-',t,depoint,'*-')
% % hold on
% % plot(tao,newpoint,'o')
% % hold off

% % error=sum((ideal-depoint).^2)
% % error=sum((ideal-depoint1).^2)
% plot(t,ideal,'o-',t,depoint,'*-',t,depoint1,'^-');
