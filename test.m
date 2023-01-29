clear
close all
clc
%%
Sys=TA_System('Sys',1e5,100,"Helium","Water",1,[100,0.5,300]);
M=TA_Minor_Loss('M',Sys,0.5,1000);
Sys.runSystem(1);
Sys.Pressure(end)