clear
close all
clc
%%%%
inert="Air";
reactive="Water";
%always two gases, one inert and one reactive
dry_switch=0;
P_m=1e5;
T_m=300;
Array=collect_properties(inert,reactive); 
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m] = Mixture_Properties(P_m,T_m,dry_switch,Array);