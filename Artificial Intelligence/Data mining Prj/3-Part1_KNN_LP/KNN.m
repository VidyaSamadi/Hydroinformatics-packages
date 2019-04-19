clc
clear
close all
format shortG
%% Insert Data

TrainData.Recharge=[433.5 273.1 201.8 235.5 281.8 275 277.9 324.6 284 224.7 113.5]';
TrainData.Discharge=[27 31.4 33.8 38 33.8 35 33.8 37 33.8 38 39]'
TrainData.WaterLvl=[.04 0 1.08 4.47 3.16 3.19 2.98 3.15 2.97 5.31 7.42]';
X(:,1)=TrainData.Recharge ;
X(:,2)=TrainData.Discharge ;
Y(:,1)=TrainData.WaterLvl ;

%% KNN Parameters

MaxK=5;
MaxRun=5;

%% Main Prog

KNNmodel=cell(MaxK,MaxRun);
cvmodel=cell(MaxK,MaxRun);
TotalResubLoss=zeros(MaxK,MaxRun);
TotalKFoldLoss=zeros(MaxK,MaxRun);

for k=1:MaxK
    for run=1:MaxRun
    KNNmodel{k,run}=ClassificationKNN.fit(X,Y,'NumNeighbors',k);
    cvmodel{k,run}=crossval(KNNmodel{k,run});
    TotalResubLoss(k,run)=resubLoss(KNNmodel{k,run});
    TotalKFoldLoss(k,run)=kfoldLoss(cvmodel{k,run});
    end
end
clc

ResubLoss=mean(TotalResubLoss,2);
KFoldLoss=mean(TotalKFoldLoss,2);

%% Select Best K

[BestValue,BestK]=min(KFoldLoss);

disp(['Best Value = ' num2str(BestValue) ])
disp(['Best K = ' num2str(BestK) ])


%%  Results

figure;
plot(ResubLoss,'r','LineWidth',2);
hold on;
plot(KFoldLoss,'b','LineWidth',2);
legend('Resub Loss','k-Fold Loss');
xlabel('Number of Neighbors - k');
ylabel('Loss');
axis([1 5 0 1.2])

%% Test

disp('Prediction of first WaterLevel with Recharge and Discharge = [375.82 29.32] : ')
disp('(Measured Value= 1.56)')
PredictedValue=KNNmodel{BestK,1}.predict([375.82 29.32])

disp('Prediction of first WaterLevel with Recharge and Discharge = [163.42 33.82] : ')
disp('(Measured Value= 2.13)')
PredictedValue=KNNmodel{BestK,1}.predict([163.42 33.82])

disp('Prediction of first WaterLevel with Recharge and Discharge = [124.16 35.12] : ')
disp('(Measured Value= 6.01)')
PredictedValue=KNNmodel{BestK,1}.predict([124.16 35.12])


