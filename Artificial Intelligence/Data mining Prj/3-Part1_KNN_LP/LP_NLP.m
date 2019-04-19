clc;
clear;
close all;

TrainData.Recharge=[433.5 273.1 201.8 235.5 281.8 275 277.9 324.6 284 224.7 113.5]';
TrainData.Discharge=[27 31.4 33.8 38 33.8 35 33.8 37 33.8 38 39]';
TrainData.WaterLvl=[.04 0 1.08 4.47 3.16 3.19 2.98 3.15 2.97 5.31 7.42]';
x(:,1)=TrainData.Recharge ;
x(:,2)=TrainData.Discharge ;
y(:,1)=TrainData.WaterLvl ;

    p=fit([x(:,1),x(:,2)],y,'poly22')
    yhat=feval(p,[x(:,1),x(:,2)]);
    
    R2=GetR2(y,yhat);
    disp(['R^2 = ' num2str(R2)]); 
    
disp('Prediction of first WaterLevel with Recharge and Discharge = [375.82 29.32] : ')
disp('(Measured Value= 1.56)')
PredictedValue=feval(p,[375.82,29.32])

disp('Prediction of first WaterLevel with Recharge and Discharge = [163.42 33.82] : ')
disp('(Measured Value= 2.13)')
PredictedValue=feval(p,[163.42,33.82])

disp('Prediction of first WaterLevel with Recharge and Discharge = [124.16 35.12] : ')
disp('(Measured Value= 6.01)')
PredictedValue=feval(p,[124.16 35.12])
    
    
% figure;
% plot(x,y,'o','MarkerSize',7);
% hold on;
% plot(x,yhat,'r','LineWidth',2);
% legend('Actual Data','Model Output');
% title(['R^2 = ' num2str(R2)]); 
