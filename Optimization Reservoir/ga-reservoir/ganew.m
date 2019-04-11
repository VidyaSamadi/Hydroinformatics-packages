

tic  

%------初始格式化--------------------------------------------------
clear all;
clc;
ps=200; MaxDT=100; D=24;  %population , tekrar , parametr
%-------------------初始化种群的个体
zmin=[970,970,970,970,952,952,952,952,952,952,970,970,888,888,888,888,888,893,888,888,888,888,888,888];
zmax=[977,977,977,980,980,977,977,957,977,970,977,977,898,898,898,898,898,898,893,892,892,898,898,898];

for i=1:ps              
    for j=1:24
    x(i,j)=zmin(j)+(zmax(j)-zmin(j))*rand;  
    end
end

for i=1:ps               %计算各个粒子的适应度，也即梯级总发电量
    f(i)=gafx(x(i,:));             % 第i个个体对应的梯级发电量
end

%------------------------------------循环计算开始第24行到第175行
    for t=1:MaxDT        %MaxDT最大迭代次数100
  
        %%---%%----------适应度求和，以便轮盘赌，初始水位为x(i,:)
            %----求累计发电量
            s(1)=f(1);
            for i=1:ps-1
                s(i+1)=s(i)+f(i+1);   %累计发电量
            end   
            %----开始轮盘赌
            for i=1:ps
                p(i)=f(i)/s(ps);   %概率
            end
            %----求累计概率
            d(1)=p(1);
            for i=1:ps-1
                d(i+1)=d(i)+p(i+1);  %累计概率
            end
            %----开始选择
            for i=1:ps            %%%%%%%%%%%%%%
                a=rand(1,ps);   %产生1行200列的一个随机数组（每个数都在0--1之间）
                ff(i)=0;        %清空并赋以初始值零
                for k=1:ps
                    if a(i)<=d(k)        %见ga算法--深入讲解第7页（此步可以改进）
                        ff(i)=f(k);      %将离得最近的值更换 ff(i)代替以前的f(i)     
                        xx(i,:)=x(k,:);  %新水位值xx(i,:)代替以前的x(i,:)
                        k=ps;            %等价于break命令，跳出循环
                    end
                end
            end
            %%--------------------------至此轮盘赌选择结束，新水位为xx(i,:)
            
       %%-----%% ---------------------下面为交叉，此时水位为xx(i,:)        
            %----交叉前准备，由大到小排序
            for i=1:ps-1    
                for j=i+1:ps          %冒泡排序法
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
            %----交叉前保存前10% 最大值    
            for i=1:0.1*ps 
                xxxx(i,:)=xx(i,:);   %xxxx(i,:)用于保存前10%的xx(i,:)
                ffff(i)=ff(i);
            end
     
            %----求群体平均值和最大值     
            favg=mean(ff);     %求平均值
            fmax=max(ff);
            %----打乱排好的顺序     
            sx=randperm(ps);   %%%%打乱产生新的数列 产生一个总数为200的数列，但是是随机产生的，杂乱数列randperm(3)可能产生 3 1 2 或者 2 3 1
            for i=1:ps    
                j=sx(i);                %新数列的下标
                ffa(i)=ff(j);           %新电量值，也即打乱原有的ff(j)的顺序，产生新的序列ffa(i)  
                xxa(i,:)=xx(j,:);       %新水位值xxa(i,:)
            end
            %--%---开始交叉（两个粒子交叉），新水位值xxa(i,:)
            for i=1:2:ps              %%两两交叉，这样可以避免混乱，ps一分为二
                if ffa(i)-ffa(i+1)>0
                    fa=ffa(i);
                else
                    fa=ffa(i+1);
                end
                if fa<favg        %%%适应度低于平均值对应0.9，适应度高对应小于0.9的pc
                    pc=0.9;
                else
                    pc=0.9-0.3*(fa-favg)/(fmax-favg);
                end
                m=rand(1,ps);  %产生1行200列的0--1的随机数字 
                if m(i)<pc          %%判断条件？？？？
                    b=rand*D+1;     %确定随机交叉位置%--D为24
                    weizhi=floor(b);   %%floor即为取整函数，便于计算
                    apc1=xxa(i,weizhi);
                    apc2=xxa(i+1,weizhi);
                    a=rand;
                    xxa(i,weizhi)=a*apc1+(1-a)*apc2;
                    xxa(i+1,weizhi)=a*apc2+(1-a)*apc1;
                end
             end              %交叉结束
             %----求交叉后的适应度值
             for i=1:ps
                ffa(i)=gafx(xxa(i,:));    %发电量
             end
             %----对交叉后的个体进行适应度排序
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
             %----将交叉后的后10%(20个)替换为先前交叉前保存的优秀个体 （20个）
            for i=0.9*ps+1:ps    
                ffa(i)=ffff(i-0.9*ps);
                xxa(i,:)=xxxx(i-0.9*ps,:);
            end
            %%---至此，交叉操作完成
            
        %%--%%---------------------开始变异
           %----将交叉个体转换成变异个体，做前期准备 
            for i=1:ps   
                 fc(i)=ffa(i);
                 xc(i,:)=xxa(i,:);          %%新的水位值xc(i,:)；；；为了储存现有数组以便备份
            end
            %------------------------与变异后的个体做比较
            for i=1:ps 
                 ffc(i)=fc(i);
                 xxc(i,:)=xc(i,:);         %新的水位值xxc(i,:)
            end
            %-----开始变异  （单个粒子变异）---新的水位值xxc(i,:)
            favg=mean(ffc);  
            fmin=min(ffc);
            %----前80%一种操作
            for i=1:ps*0.8     
               if fc(i)<=favg    %选择函数
                  pm=0.1;
               else
                  pm=0.1-0.099*(fc(i)-favg)/(favg-fmin);   %%自己发挥的
               end
               b=rand;
               if b<pm
                  j=floor(24*rand)+1;
                  xxc(i,j)=zmin(j)+(zmax(j)-zmin(j))*rand;
               end
            end
            %----后80%一种操作---梯度投影法
            for i=(ps*0.8+1):ps            
                  [mem,pos]=max(ffc);      %pos为最大的ffc对应的下标值   
                  temp=xxc(pos,:);         %temp为ffc的最大值
                  v=temp-xc(i,:);  
                  xc(i,:)=xc(i,:)+rand*v;  %--rand*v为单位梯度，从xc(i,:)指向temp
                  fc(i)=gafx(xc(i,:));
                        if  fc(i)>ffc(i)
                            ffc(i)=fc(i);
                            xxc(i,:)=xc(i,:);
                        end
            end
            %%----至此变异操作结束
            
       %%--%%-----------为下一代做数据替换
       for i=1:ps
           f(i)=ffc(i);
           x(i,:)=xxc(i,:);  %最终的水位值xxc(i,:)
       end
           
    end
  
    %-------------------------------------至此一个遗传循环结束
    %--------------------------------下一段为分类输出重新计算
    [mem,pos]=max(ffc);
    temp=xxc(pos,:) ;  %最大电量对应的水位
     x=temp            %最优水位值
     y=mem             %最大发电量
    toc     %-----此时记一次时间
    
    