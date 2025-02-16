clear all;
close all;
clc;

% x = 1;
% for i = i = 5:7:62
% datax = readmatrix('SiN_SiO2__Channel2.25.xlsx');
% % eff(1,x) = data(12,i+1); 	  	% effective mode area [m^2]
% % rows = 1:81;  % Select row indices
% % cols = i ;  % Select column indices
% % n_eff(:,x) = data(rows, cols); % Extract specific rows and columns:
% Wc(1,x) = datax(1,i-6);
% Wc1(1,x) = datax(1,i-5);
% Hc(1,x) = datax(1,i-4);
% Hc1(1,x) = datax(1,i-3);
% x = x+1;
% 
% end

x = 1;

for i = 5:7:62
    datax = readmatrix('SiN_SiO2__Channel2.25.xlsx'); % Read entire Excel file
    Wc(1,x) = datax(1,i-4);
    Hc(1,x) = datax(1,i-3);
    eff(1,x) = datax(12,i+1);
    x = x+1;
end