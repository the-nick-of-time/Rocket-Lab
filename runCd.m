clear
close all
clc
path(path, 'C:\Users\Nicholas\Documents\MATLAB\MATLAB-customs')

%files = {'./TUNNEL/Group5_30Mar.txt', './TUNNEL/Group8_Fin1.txt', './TUNNEL/Group8_Fin2.txt', './TUNNEL/Group14_test1.txt', './TUNNEL/Group17_test1.txt'};
files = {'./TUNNEL/Group17_test1.txt'};
data = windtunneldata(files);

[cd, stdcd] = calculateCd(data, .008659)