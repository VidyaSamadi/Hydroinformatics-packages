function[varargout]=AWBM(c1,c2,c3,a1,a2,kb,bfi)

%% Global variables
global n et prt

%% pre allocated matrix size
qt = zeros(n,1);    bst = zeros(n,1);   s1t = zeros(n,1);
exs = zeros
exs(1)=0;
qbt(1)=0;
qt(1)=0;
exb(1)=0;

a3=1-a1-a2;
for i = 2:n
%%	storage 1
    pet1=((prt(i)+ s1t(i-1))/c1) * et(i);
    if pet1 > et(i) 
       pet1=et(i); 
    end

	s1t(i) = s1t(i-1) + (prt(i) - pet1);
        if s1t(i)<=0 
            s1t(i)=0;
        
        end
    
    if s1t(i) < c1 
       ex1(i) = 0;
    else
       ex1(i) = s1t(i) - c1;
       s1t(i) = s1t(i) - ex1(i);
    end
    
%%	storage 2
    pet2=((prt(i)+ s2t(i-1))/c2) * et(i);
    if pet2 > et(i)
       pet2=et(i);
    end
		
    s2t(i) = s2t(i-1) + (prt(i) - pet2);
    if s2t(i)<= 0
       s2t(i) = 0;
    end
    
    if s2t(i) < c2
       ex2(i) = 0;
    else
       ex2(i) = s2t(i) - c2;
       s2t(i) = s2t(i) - ex2(i);
    end
    
%%	storage 3
    pet3=((prt(i)+ s3t(i-1))/c3) * et(i);
     if pet3 > et(i) 
        pet3=et(i);
     end

    s3t(i) = s3t(i-1) + (prt(i) - pet3);
    if s3t(i) <= 0
       s3t(i) = 0;
    end
    
    if s3t(i) < c3
       ex3(i) = 0;
    else
       ex3(i) = s3t(i) - c3;
       s3t(i) = s3t(i) - ex3(i);
    end
					
%%	excess runoff from storage
		ex(i) = (ex1(i)*a1) + (ex2(i)*a2) + (ex3(i)*a3);
%%	surface runoff
		exs(i) = (1 - bfi) * ex(i);
%%	baseflow
		exb(i) = bfi * ex(i);
        
	    bst(i)=bst(i-1) + exb(i);
            if bst(i)<0.001
                bst(i)=0;
            
            end
        
		qbt(i) = (1-kb)*bst(i);
        
		bst(i) = bst(i) - qbt(i);
            if bst(i)<0
                bst(i) = 0;
        
            end


%%	total flow
	qt(i) = exs(i) + qbt(i);
end

varargout{1} = qt;
varargout{2} = exs;
varargout{3} = qbt;
