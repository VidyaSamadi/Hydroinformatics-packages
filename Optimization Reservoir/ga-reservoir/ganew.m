

tic  

%--------------------------------------------------------
clear all;
clc;
ps=200; MaxDT=100; D=24;  %population , tekrar , parametr
%-------------------
zmin=[970,970,970,970,952,952,952,952,952,952,970,970,888,888,888,888,888,893,888,888,888,888,888,888];
zmax=[977,977,977,980,980,977,977,957,977,970,977,977,898,898,898,898,898,898,893,892,892,898,898,898];

for i=1:ps              
    for j=1:24
    x(i,j)=zmin(j)+(zmax(j)-zmin(j))*rand;  
    end
end

for i=1:ps               %
    f(i)=gafx(x(i,:));             % 
end

%------------------------------------24,175
    for t=1:MaxDT        %MaxDT 100
  
        %%---%%----，，x(i,:)
            %----
            s(1)=f(1);
            for i=1:ps-1
                s(i+1)=s(i)+f(i+1);   
            end   
            %---
            for i=1:ps
                p(i)=f(i)/s(ps);   %
            end
            %----
            d(1)=p(1);
            for i=1:ps-1
                d(i+1)=d(i)+p(i+1);  %
            end
            %----
            for i=1:ps            %%%%%%%%%%%%%%
                a=rand(1,ps);   2000--1
                ff(i)=0;        %
                for k=1:ps
                    if a(i)<=d(k)        %
                        ff(i)=f(k);      % ff(i)f(i)     
                        xx(i,:)=x(k,:);  %xx(i,:)x(i,:)
                        k=ps;            %
                    end
                end
            end
            %%--------------------------xx(i,:)
            
       %%-----%% ---------------------xx(i,:)        
            %----
            for i=1:ps-1    
                for j=i+1:ps          %
                   if ff(i)<=ff(j)
                        b=ff(i);
                        c=xx(i,:);
                        ff(i)=ff(j);
                        xx(i,:)=xx(j,:);
                        ff(j)=b;
                        xx(j,:)=c;  
                    end
                end
            end    
            %---10%     
            for i=1:0.1*ps 
                xxxx(i,:)=xx(i,:);   %xxxx(i,:)10%的xx(i,:)
                ffff(i)=ff(i);
            end
     
            %----   
            favg=mean(ff);     %
            fmax=max(ff);
            %----   
            sx=randperm(ps);   %%%%200，，randperm(3) 3 1 2  2 3 1
            for i=1:ps    
                j=sx(i);                %
                ffa(i)=ff(j);           %ff(j)ffa(i)  
                xxa(i,:)=xx(j,:);       %xxa(i,:)
            end
            %--%---，xxa(i,:)
            for i=1:2:ps              %%，ps
                if ffa(i)-ffa(i+1)>0
                    fa=ffa(i);
                else
                    fa=ffa(i+1);
                end
                if fa<favg        %%%0.9，0.9pc
                    pc=0.9;
                else
                    pc=0.9-0.3*(fa-favg)/(fmax-favg);
                end
                m=rand(1,ps);  %12000--1
                if m(i)<pc          %%？？？？
                    b=rand*D+1;     %%--D24
                    weizhi=floor(b);   %%floor
                    apc1=xxa(i,weizhi);
                    apc2=xxa(i+1,weizhi);
                    a=rand;
                    xxa(i,weizhi)=a*apc1+(1-a)*apc2;
                    xxa(i+1,weizhi)=a*apc2+(1-a)*apc1;
                end
             end              %
             %----
             for i=1:ps
                ffa(i)=gafx(xxa(i,:));    %
             end
             %----
             for i=1:ps-1  
                for j=i+1:ps
                    if ffa(i)<ffa(j)
                        b=ffa(i);
                        c=xxa(i,:);
                        ffa(i)=ffa(j);
                        xxa(i)=xxa(j);
                        ffa(j)=b;
                        xxa(j,:)=c;
                    end
                end
             end 
             %----10%(20) （20）
            for i=0.9*ps+1:ps    
                ffa(i)=ffff(i-0.9*ps);
                xxa(i,:)=xxxx(i-0.9*ps,:);
            end
            %%--
            
        %%--%%-------------------
           %---
            for i=1:ps   
                 fc(i)=ffa(i);
                 xc(i,:)=xxa(i,:);          %%xc(i,:)；；；
            end
            %------------------------
            for i=1:ps 
                 ffc(i)=fc(i);
                 xxc(i,:)=xc(i,:);         %xxc(i,:)
            end
            %--------xxc(i,:)
            favg=mean(ffc);  
            fmin=min(ffc);
            %----80%一
            for i=1:ps*0.8     
               if fc(i)<=favg    %
                  pm=0.1;
               else
                  pm=0.1-0.099*(fc(i)-favg)/(favg-fmin);   %%
               end
               b=rand;
               if b<pm
                  j=floor(24*rand)+1;
                  xxc(i,j)=zmin(j)+(zmax(j)-zmin(j))*rand;
               end
            end
            %----80%一
            for i=(ps*0.8+1):ps            
                  [mem,pos]=max(ffc);      %pos  
                  temp=xxc(pos,:);         %tempffc
                  v=temp-xc(i,:);  
                  xc(i,:)=xc(i,:)+rand*v;  %--rand*v，从xc(i,:)temp
                  fc(i)=gafx(xc(i,:));
                        if  fc(i)>ffc(i)
                            ffc(i)=fc(i);
                            xxc(i,:)=xc(i,:);
                        end
            end
            %%----
            
       %%--%%----------
       for i=1:ps
           f(i)=ffc(i);
           x(i,:)=xxc(i,:);  %xxc(i,:)
       end
           
    end
  
    %-------------------------------------
    %--------------------------------
    [mem,pos]=max(ffc);
    temp=xxc(pos,:) ;  %
     x=temp            %
     y=mem             %
    toc     %-----
    
    
