clear
clc
%% load data
Data
ReadFromCSV
Dimensionless
%%
M = [eye(5), CB.CB_MK.Max;
     CB.CB_MK.Max', CB.CB_MK.Mxx];
K = [diag(CB.CB_MK.Kaa), zeros(5,12);
     zeros(12,5), CB.CB_MK.Kxx];
F = [CB.CB_F.Fa;
     CB.CB_F.Fx];