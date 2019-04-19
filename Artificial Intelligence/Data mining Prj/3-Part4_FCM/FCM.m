clc
clear
close all

x=[17.89 24.26 38.91 18.67 20.67 23.11 35.55 32.14 ...
   22.53 18.58 16.56 14.11 13.42 13.06 12.62 16.69 17.63 26.43];

disp('Number of Clusters = 2 ' )


[center_2,U_2,Objfunc_2]=fcm(x',2);
[center_3,U_3,Objfunc_3]=fcm(x',3);

disp('Center of clusters:')
center_2

disp('Membership Function Matrix:' )
U_2

disp('Number of Clusters = 3 ' )

disp('Center of clusters:')
center_3

disp('Membership Function Matrix:' )
U_3

