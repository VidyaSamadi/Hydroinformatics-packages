clc;
clear;
close all;

n=100;

x=[10 30 50 70 110 140 170 200]';
X=[x ones(length(x),1)];
Y=[105 188 200 208 233 244 253 260]';



theta=pinv(X)*Y;
disp(['y=' num2str(theta(1)) '*x+(' num2str(theta(2)) ')'])
xmin=min(x);
xmax=max(x);
xx=linspace(xmin,xmax,1000);
yy=theta(1)*xx+theta(2);
xn=[100 125 170 60]
Yhat(1)=polyval(theta,xn(1));
Yhat(2)=polyval(theta,xn(2));
Yhat(3)=polyval(theta,xn(3));
Yhat(4)=polyval(theta,xn(4))
figure;
plot(x,Y,'o');
hold on;
plot(xn,Yhat,'o','MarkerFaceColor','g');
hold on;
plot(xx,yy,'r','LineWidth',2);
legend('Data','Linear Model (LS)');
xlabel('x');
ylabel('y');
title('Linier Regression')


