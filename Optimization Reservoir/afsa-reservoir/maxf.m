function result = maxf(x)
%-----------------------------------
%--------��������------------------
%%%%%%%%%%%%%%%%% ÿ����ˮ������ˮ�����
qwan=[ 313 	363 	587 	493 	263 	309 	699 1110 	1134 	808 	486 	325 ];
%qwan=[402,352,702,457,227,432,418,1040,795,544,373,405]; 
%qwan=[598.21,699.12,564.58,667.52,195.60,195.60,195.60,946.89,3235.13,1505.35,951.89,730.58]; %��ˮ��
%qwan=[167.16,176.81,378.46,317.67,286.8,286.8,1061.56,1312.18,918.52,768.25,500.99,250.13];%ƽˮ��
%qwan=[332.3,295.18,664.44,365.42,26.15,45.32,122.58,950.26,622.13,141.65,236.98,334.18];%��ˮ��
%%%%%%%%%%%%%%% ÿ����ˮ����ˮ�����
qwl=[0,0,0,0,0,0,0,0,0,0,0,0];  %������ˮ����������w�������կ��l��������

%%%%%%%%%%%%ÿ��Сʱ��
t(1)=31*24; t(3)= t(1); t(5)= t(1);
t(5)= t(1); t(7)= t(1); t(8)= t(1);
t(10)= t(1); t(12)= t(1);
t(2)=28*24;
t(4)=30*24; t(6)= t(4);t(9)= t(4);t(11)= t(4);

%-------------------------------------------------

     for j=1:12        %%%%%%% ˮλ�������
        if 952<=x(j)<955
            vw(j)=(4.508+(x(j)-952)/(955-952)*(4.548-4.508))*10^8;  %vw(j)����ˮλ��Ӧ�Ŀ���
        elseif 955<=x(j)<957
            vw(j)=(4.548+(x(j)-955)/(957-955)*(4.739-4.548))*10^8;
        elseif 957<=x(j)<960
            vw(j)=(4.739+(x(j)-957)/(960-957)*(4.936-4.739))*10^8;
        elseif 960<=x(j)<965
            vw(j)=(4.936+(x(j)-960)/(965-960)*(5.546-4.936))*10^8;
        elseif 965<=x(j)<970
            vw(j)=(5.546+(x(j)-965)/(970-965)*(7.756-6.563))*10^8;
        elseif 975<=x(j)<977.5
            vw(j)=(7.756+(x(j)-975)/(977.5-975)*(8.35-7.756))*10^8;
        elseif 977.5<=x(j)<980
            vw(j)=(8.35+(x(j)-977.5)/(980-977.5)*(8.962-8.35))*10^8;
        end
     end
 
             
     for j=1:11    %%%%%%%%%%%%%%% ��й������ˮλ
         qxw(j)=qwan(j)-(vw(j+1)-vw(j))/(t(j)*3600); % qxw(j)Ϊ��������
     end
         qxw(12)=qwan(12)-(vw(1)-vw(12))/(t(12)*3600);  %vw(13)=vw(1)
     %%%%%%%%%%%%%%%%%%%%%%%%%% ��й����ˮλ��ϵ
     for j=1:12
         if 0<=qxw(j)<85
             hxw(j)=898+(qxw(j)-0)/(85-0)*(899-898);  %����������Ӧ��ˮλ
         elseif 85<=qxw(j)<174
             hxw(j)=899+(qxw(j)-85)/(174-85)*(900-899);
         elseif 174<=qxw(i)<286
             hxw(j)=900+(qxw(j)-174)/(286-174)*(901-900);
         elseif 286<=qxw(j)<571
             hxw(j)=901+(qxw(j)-286)/(571-2860)*(902-901);
         elseif 571<=qxw(j)<878
             hxw(j)=902+(qxw(j)-571)/(878-571)*(903-902);
         elseif 878<=qxw(j)<1320
             hxw(j)=903+(qxw(j)-878)/(1320-878)*(904-903);
         elseif 1320<=qxw(j)<1860
             hxw(j)=904+(qxw(j)-1320)/(1860-1320)*(905-904);
         elseif 1860<=qxw(j)<2480
             hxw(j)=905+(qxw(j)-1860)/(2480-1860)*(906-905);
         end
     end
      %%%%%%%%%%%%%%%%%%%%%%%%%%% ������ˮλ��
      for j=1:11
          hw(j)=(x(j)+x(j+1))/2-hxw(j); %����ƽ��ˮλ-����ˮλ  %%%%%% �˴��޸�
      end
          hw(12)=(x(12)+x(1))/2-hxw(12);
      %%%%%%%%%%%%%%%%%%%%%%% ����������ȡA=8.0
      for j=1:12
          nw(j)=8.0*qxw(j)*hw(j);  %����ϵ��ȡΪ8.0�����ԸĽ�����ÿ���³���ֵ 
          ew(j)=8.0*qxw(j)*hw(j)*t(j); %��ÿ���µ���ֵ
      end
      sw=ew(1); 
      for j=2:12
          sw=sw+ew(j);  % �����ۼ� ���
      end
    
  
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%�ݼ���վ�������귢����
        result=sw;        %�����ܷ������� 
