function [NTrain,NTest] = Normalize(Train,Test)
% Normalize
a=numel(Train(1,:));
for i=1:a
    NTrain(:,i)=(Train(:,i)-min(Train(:,i)))/(max(Train(:,i))-min(Train(:,i)));
    NTest(:,i)=(Test(:,i)-min(Train(:,i)))/(max(Train(:,i))-min(Train(:,i)));
end

