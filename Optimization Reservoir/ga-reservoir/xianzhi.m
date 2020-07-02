function [XX] =xianzhi(XXnext)  %
format short    
    a=[977,977,977,980,980,977,977,957,977,970,977,977];   %%[898,898,898,898,898,898,893,892,892,898,898,898];   %
    b=[970,970,970,970,952,952,952,952,952,952,970,970];   %%[888,888,888,888,888,893,888,888,888,888,888,888];   %
    for i=1:12   
        if(XXnext(i)>a(i))  %XXnext1(1)qq1，XXnext1(2)qq2
        XXnext(i)=a(i);   
        end
        if(XXnext(i)<b(i))
        XXnext(i)=b(i);   
        end
    end
    XX=XXnext;
