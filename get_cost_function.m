function Fvalue=get_cost_function(kappa)
%% option 1
M_array=collect_properties2('DimethyEther','Water',kappa);
calualtedv=Mixture_Properties(pm,T,1,M_array);
%%
measuredv=49324;
Fvalue=max(abs(calualtedv-measuredv));
end