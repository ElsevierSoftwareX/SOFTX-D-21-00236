function [k,cp,rho] = solid_properties(solid,T)
%{
this function calculates solid properties based
on the temperature (k) and solid type
all results are in SI units. all results are based
on the data in DeltaEC Users Guide
Version 6.4b2.7 unless said otherwise
% results were verified against deltaEC
%}

%solid_properties("Copper",300) will calculate the properties of Copper at
%300k

%solid_properties({"Copper",0.4},300) will calculate the properties of Copper at
%300k, and perform a weighted average with the properties of an ideal
%solid. this allows to change the material smoothly
%

%% if solid is with a ratio
if length(solid)==2&&str2double(solid(2))<=1&&str2double(solid(2))>=0
    ratio=str2double(solid(2));
    [k1,cp1,rho1]=solid_properties(solid(1),T);
    k=k1*ratio+1e8*(1-ratio);
    cp=cp1*ratio+1e8*(1-ratio);
    rho=rho1*ratio+1e8*(1-ratio);
    return
elseif length(solid)~=1
    error('wrong input format')
end

switch solid
    case 'Ideal'
        k=1e8;
        rho=1e8;
        cp=1e8;
    case 'Copper' 
        k=398-0.0567*(T-300);
        rho=9000;
        cp=420;
    case 'Nickel' 
        if T<631
            k=63.8+0.08066*(631-T);
        else
            k=63.8+0.02156*(T-631);
        end
        rho=8700;
        cp=530;   
    case 'Stainless'
        rho=8274.55-1055.23*exp(-((T-273.15-2171.05)/2058.08)^2);
        k=(266800*T^(-5.2)+0.21416*T^(-1.6))^-0.25;
        cp=(1.7054e-6*T^(-0.88962)+23324/T^6)^(-1/3)+15/T;
        
    case 'Molybdenum'
        rho=10868.6-2637.52*exp(-((T-273.15-11383.7)/9701.36)^2);
        k=(33.9616-0.00947953*(T-273.15)-4.12809e-08*(T-273.15)^2)*4.186;
        cp=253.791+0.0583812*(T-273.15)-2.73919e-06*(T-273.15)^2;
    case 'Tungsten'
        rho=19254*(1-3*(-8.69e-5+3.83e-6*(T-273.15)+7.92e-10*(T-273.15)^2));
        k=135.5+1.05e4/T-0.023*T;
        cp=0.13576e3*(1-4805/T^2)+0.0091159*T+2.31341e-9*T^3;
    case 'Mylar'
        rho=1400-0.175*T;
        k=0.11+1.7e-4*T;
        cp=3.7*T;
    case 'Kapton'
        rho=1445-0.085*T;
        k=0.2*(1-exp(-T/100));
        cp=3.64*T;
    case 'Celcor'
        rho=2510;
        k=2.5;
        cp=262.5+1.864*T-0.001011*T^2;
    case 'Ceramic' %ceramic stack
        %data is from MatWeb material property data 2021
        k= 2.5;
        rho=2300;
        cp= 900;

    otherwise
        error('no data for specified solid. check for spelling and case')
end
if isnan(k)||isnan(rho)||isnan(cp)
    error 'nan'
end
end

