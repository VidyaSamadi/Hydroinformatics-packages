%% Calculate the weighted covariance matrix
function xy = covweight(x,w)

[m,n] = size(x);

W=zeros(size(x,1),size(x,2))*NaN;
for repeat=1 : size(x,2)
W(:,repeat)=w';
end

center1=W.*x;

 
  xc = bsxfun(@minus,x,sum(center1));  
  
Wx=zeros(size(x,1),size(x,2))*NaN;
for repeat=1 : size(x,2)
Wx(:,repeat)=sqrt(w');
end

xy2=Wx.*xc;

xy=xy2'*xy2;
  
    
    
 

