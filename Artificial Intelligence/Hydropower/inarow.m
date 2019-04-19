function [ Y ] = inarow( X )
y=[]
for i=1:size(X,1)
    if i==1
        Y=X(i,:);
    else
    Y= [Y , X(i,:)];
    end
end

