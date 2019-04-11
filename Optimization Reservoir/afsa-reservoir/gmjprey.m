function [XXnext] = gmjprey(XX,Try_number,Visual,Step)  %觅食行为 
pp=0;
for j=1:Try_number
    XXj=XX+rand*Step*Visual;
    if(maxf(XX)<maxf(XXj))
        XXnext=XX+rand*Step*(XXj-XX)/norm(XXj-XX);  %二范数
        pp=1;
        break
    end
end
if(~pp)  %即为pp不等于0的时候
   XXnext=XX+rand*Step;
end