function [ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m] = Mixture_Properties(pm,T,dry_switch,M_Array,varargin)
%% function description
%This function calculates the properties of a gas mixture based on
% temperature, pressure, inert and reactive gas. the mixture is considered
%to be in equilibrium with a liquid layer of the reactive component, or
%else dry (if dry_switch=1)

%the function is valid in the pressure range between 1 and 10 bars. it
%returns values for pressures between 0.9 and 20 bars, but these should not
%be considered accurate if out of range.

% The function is based on the book "the properties of gases and liquids" by
% Bruce Poling, fifth edition (2001) McGraw-Hill. all sections referenced
% in this code are from the same book.

%% input parameters:
%pm         - mixture presure [pa]
%T          - mixture temperature  [k]

%dry_switch - a parameter indicating if the system includes phase change
%should be 1 for the dry case and 0 for the wet case.
%if a value between 0 and 1 is input for dry_switch,
%the function will calculate the dry and wet case seperatly
%and average the results based on the value of dryswitch.
%this option has no physical meaning, it is used to transfer
%between dry and wet modes without diverging.

%M_Array        a cell array containing relevant data from the property
% tables. the data should be collected from the proerty
% tables using the  function "collect_properties".

%if the mixture inludes several reactives, to be decided between based on
%condtitions, the last input, to be input in varargin, should be a
%function handle for the decision function f(pm,T) returning two indexes
%based on pressure and temperature
%for example, if the properties were collected using:
%M_Array=collect_properties("Air",["Water";"Acetone"])
%the decision function can be
%f=@(pm,T)[1,(T>350)+1]
%Mixture_Properties(pm,T,dry_switch,M_Array,f)
%and the properties function will use water for temperatures above 350 and
%acetone otherwise.
%% output parameters:
% nu      - kinematic viscocity          [m^2/s]
% alpha   - thermal diffusivity          [m^2/s]
% D       - diffusion coefficient        [m^2/s]
% rho     - density                      [kg/m^3]
%Pr       - prandtl number               [-]
%Sc       - schmidt number               [-]
%gamma    - specific heat ratio          [-]
%M        - mixture molar mass           [g/mol]
%cp_mol   - molar heat capacity          [j/mol k]
%cp       - heat capacity                [j/kg k]
%k        - thermal conductivity         [W/m k]
%sp       - speed of sound               [m/s]
%lh       - Moldar latent heat           [kj/mol]
%Tref     - reactive boiling temperature [k]
%CB       - mass fraction of reactive    [-]
%lh_m     - latent heat in kj/kg
%% special cases

%if the reactive or inert is a mixture
if length(M_Array)==2
    [ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m] =...
        Double_Mixture(pm, T,dry_switch,M_Array{1,1},M_Array{1,2},M_Array{2,1},M_Array{2,2});
    return
end

% if the reactive or inert are a "decision mixture"
if length(size(M_Array))==3
    if isempty(varargin)
        error('no decision function supplied to decide between gases')
    elseif length(varargin)>1||~isa(varargin{1},'function_handle')
        error('wrong format for decision function')
    else
        f=varargin{1};
        idxs=f(pm,T);
        M_Array=M_Array(:,idxs(1),idxs(2));
    end
end

%semi-dry mixture
if dry_switch>0&&dry_switch<1
    %warning('dry switch is between 0 and 1, returning averaged value')
    [ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m] =...
        Semi_Dry_Mixture(pm,T,dry_switch,M_Array);
    return
end
%% check for errors
if pm>10e6||pm<9e4              %if the pressure is out of working range
    error ('pressure out of range')
end
if isnan(M_Array{1})&&~dry_switch %if the selected reactive component is inert
    error ('selected reactive does not change phase. CC equation cannot be used')
end
%% define parameters
MA = M_Array{2};      %molar masses                  [g/mol]
MB = M_Array{3};

miuA = M_Array{4};    %dipole moments                [debyes]
miuB = M_Array{5};

TbA = M_Array{6};     %boiling temperatures          [k]
TbB = M_Array{7};

TcA = M_Array{8};     %critical temperatures         [k]
TcB = M_Array{9};

omegaA = M_Array{10}; %accentric factors             [-]  
omegaB = M_Array{11};

kA = M_Array{12};     %special correction factor     [-]                
kB = M_Array{13};     %for highly polar substances, 
                    %hydrogen or helium
                    
aA = M_Array{14};     %coefficiencts for             [-]
aB = M_Array{15};     %caclulating heat capacity

VcA = M_Array{16};    %critical volumes              [cm^3/mol]
VcB = M_Array{17};
%% Calculate mole fraction using Clausius-Clapeyron equation
c_lh=M_Array{18};                                 %coefficients for LH
lh_m=(c_lh(1)*(T/500).^c_lh(2)+c_lh(3))*1e6;      %calculate latent heat in j/kg

c_tb=M_Array{19};                                 %polynomial for TB
Tref=(c_tb(1)*(pm/10^5).^c_tb(2)...               %calculate Tb [k]
    +c_tb(3)*pm/10^5+c_tb(4))*300;  

lh_b_m=(c_lh(1)*(Tref/500).^c_lh(2)+c_lh(3))*1e6; %TB latent heat in j/kg
L_H_B = lh_b_m *MB / (8314 * Tref);               %scaled TB latent heat

%CC equation
    CB = exp(-L_H_B * (Tref/T - 1)); 
if dry_switch
    CB=0;
elseif CB>1
    error('Mixture_properties:aboveB','working temperature is above boiling temperature for the reactive component')
end

CA = 1 - CB;

%% Molecular Diffusion Coefficient
%calculation method is from section 11-3 (page 639) in poling's book

% P in the diffusion calculation is given in [atm] 
P = pm / 101325;

if miuA == 0 && miuB == 0  % Non-polar gases
    
    %Lennard-Jones parameters- epsilon and sigma
    ekA = M_Array{20};    %epsilon/k [k]
    ekB = M_Array{21};
    sigA = M_Array{22};   %sigma [angstram]
    sigB = M_Array{23};
    
    sigAB = (sigA+sigB)/2;
    ekAB = sqrt(ekA*ekB);
    Td = T/ekAB;
    CI = CId(Td);
    
    D = 1e-4 * 0.0018583*sqrt(T^3*(1/MA+1/MB))/(P*sigAB^2*CI);
    
else       % Polar gases
    VbA = 0.285*VcA^1.048;
    VbB = 0.285*VcB^1.048;
    
    lpA = 1.94*10^3*miuA^2/(VbA*TbA);
    ekA = 1.18*(1+1.3*lpA^2)*TbA;
    sigA = ((1.585*VbA)/(1+1.3*lpA^2))^(1/3);
    
    lpB = 1.94*10^3*miuB^2/(VbB*TbB);
    ekB = 1.18*(1+1.3*lpB^2)*TbB;
    sigB = ((1.585*VbB)/(1+1.3*lpB^2))^(1/3);
    
    ekAB = sqrt(ekA*ekB);
    sigAB = sqrt(sigA*sigB);
    lpAB = sqrt(lpA*lpB);
    
    Td = T/ekAB; % Dimensionless Temperature
    CI = CId(Td)+0.19*lpAB^2/Td; % Corrected Collision integral for polar gases.
    D = 1e-4 * 0.0018583*sqrt(T^3*(1/MA+1/MB))/(P*sigAB^2*CI);
end

%% Dynamic Viscosity
%calculation method is based on section 9.4 (page 470) in poling's book
%the method is Chung's method, based on:
%Chung, T.-H., M. Ajlan, L. L. Lee, and K. E. Starling: Ind. Eng. Chem. Res. 27: 671 (1988).

sigA = 0.809*VcA^(1/3);
sigB = 0.809*VcB^(1/3);
sigAB = sqrt(sigA*sigB);
sigm = (CA^2*sigA^3 + CB^2*sigB^3 + 2*CA*CB*sigAB^3)^(1/3);

ekA = TcA/1.2593;
ekB = TcB/1.2593;
ekAB = sqrt(ekA*ekB);
ekm = (ekA*CA^2*sigA^3 + ekB*CB^2*sigB^3 + 2*ekAB*CA*CB*sigAB^3)/sigm^3;

Td = T/ekm;

MAB = (2*MA*MB)/(MA+MB);
Mm = ((sqrt(MA)*ekA*CA^2*sigA^2 + sqrt(MB)*ekB*CB^2*sigB^2 + 2*sqrt(MAB)*ekAB*CA*CB*sigAB^2)/(ekm*sigm^2))^2;

omegaAB = (omegaA + omegaB)/2;
omegam = (omegaA*CA^2*sigA^3 + omegaB*CB^2*sigB^3 + 2*omegaAB*CA*CB*sigAB^3)/sigm^3;

kAB = sqrt(abs(kA*kB));
km = CA^2*kA + CB^2*kB + 2*CA*CB*kAB;

mium = (sigm^3*ekm*((CA^2*miuA^4)/(sigA^3*ekA) + (CB^2*miuB^4)/(sigB^3*ekB) + 2*CA*CB*(miuA^2)*(miuB^2)/(sigAB^3*ekAB)))^0.25;


Tcm = 1.2593*ekm;
Vcm = (sigm/0.809)^3;
miurm = (131.3*mium)/sqrt(Vcm*Tcm);
Fcm = 1-0.275*omegam+0.059035*miurm^4+km;

mu = 40.785*Fcm*sqrt(Mm*T)/(Vcm^(2/3)*CIv(Td));

%% Thermal Conductivity
%calculation method is based on section 10.4 (page 576) in poling's book
%the method is Chung's method, based on:
%Chung, T.-H., M. Ajlan, L. L. Lee, and K. E. Starling: Ind. Eng. Chem. Res. 27: 671 (1988).

R = 8.314; %universal gas constant j/molK
%heat capacity calculation, see appendix a section c
cpA = R*(aA(1) + aA(2)*T*10^-3 + aA(3)*(T^2)*10^-5 + aA(4)*(T^3)*10^-8 + aA(5)*(T^4)*10^-11);
cvA = cpA - R;
cpB = R*(aB(1) + aB(2)*T*10^-3 + aB(3)*(T^2)*10^-5 + aB(4)*(T^3)*10^-8 + aB(5)*(T^4)*10^-11);
cvB = cpB - R;
cvm = CA*cvA + CB*cvB;

%for non polar gases, beta is calculated based on the mixture
%for a mixture
if (miuA == 0 && miuB == 0)||(dry_switch && miuA == 0) %non polar gas
    betam = 0.7862 - 0.7109*omegam + 1.3168*omegam^2;
elseif ~miuA                                           % one polar gas
    betam= M_Array{24};
    
else                                                   %two polar gases
    if ~dry_switch
        warning('cannot calculate conductivity for two polar gases')
        warning('calculation will be inaccurate')
    end
    betam= M_Array{25};
end
alpham = cvm/R-3/2;
Trm = T/Tcm;
zm = 2 + 10.5*Trm;
psi = 1 + alpham*((0.215+0.28288*alpham - 1.061*betam + 0.26665*zm)/(0.6366 + betam*zm + 1.061*betam));

k = 3.75*(mu*10^-7)*R*psi/(Mm/1000);

%% Heat Capacities

cp_mol = CA * cpA + CB * cpB;
cv_mol = CA * cvA + CB * cvB;

M = CA * MA + CB * MB;
cp = 1e3 * cp_mol / M;
gamma = cp_mol / cv_mol;

%% Kinematic Viscosity, Thermal Diffusivity, Pr, Sc,Sound speed, lh
rho = pm * M / (1e3 * R * T);
nu = 1e-7 * mu / rho;
alpha = k / (rho * cp);
Pr = nu / alpha;
Sc = nu / D;
sp=(gamma*pm/rho)^0.5;
lh=lh_m*MB/1000;
%% check for nans or negatives
if max(isnan([ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,...
        cp,k,sp,lh_m,MB,Tref,CB]))||max([ nu, alpha, D, rho, Pr, Sc, ...
        gamma, M, cp_mol,cp,k,sp,lh_m,MB,Tref,CB])<0
    error ('nan or negative number in parameters, something went wrong')
end
end





%% internal functions
function CI = CId(Td)
%Calculates the Collision Integral of Diffusion

A = 1.06036; %Coefficients
B = 0.15610;
C = 0.19300;
D = 0.47635;
E = 1.03587;
F = 1.52996;
G = 1.76474;
H = 3.89411;
CI = A/(Td^B)+C/exp(D*Td)+E/exp(F*Td)+G/exp(H*Td); %Collision Integral fit
end
function CI = CIv(Td)
%Calculates the Collision Integral of Viscocity and Thermal Conductivity
A = 1.16145; %Coefficients
B = 0.14874;
C = 0.52487;
D = 0.77320;
E = 2.16178;
F = 2.43787;
CI = A/(Td^B)+C/exp(D*Td)+E/exp(F*Td); %Collision Integral fit
end
function [ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m] = Double_Mixture(pm, T,dry_switch,Array1,Array2,ratio1,ratio2)
%this function calculates the properties for a "mixture" of two
%mixtures. essentially averaging them. is used to transfer between two gas
%mixtures. doesnt have a lot of physical meaning
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m]=...
    Mixture_Properties(pm, T,dry_switch,Array1);
res1=[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m];
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m]=...
    Mixture_Properties(pm, T,dry_switch,Array2);
res2=[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m];
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m]=...
    deal_array(res1*ratio1+res2*ratio2);

end
function [ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m] = Semi_Dry_Mixture(pm,T,dry_switch,Array)
%this function calculates the properties for a "semi-dry mixture"
%essentially averaging wet and dry them according to the value of
%dryswitch. used to transfer between dry and wet mode.
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m]=...
    Mixture_Properties(pm, T,1,Array);
res_wet=[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m];
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m]=...
    Mixture_Properties(pm, T,0,Array);
res_dry=[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m];
[ nu, alpha, D, rho, Pr, Sc, gamma, M, cp_mol,cp,k,sp,lh,Tref,CB,lh_m]=...
    deal_array(res_wet*dry_switch+res_dry*(1-dry_switch));

end
function varargout = deal_array(arr)
%assisting function that allows to deal an array into several named
%parameters.
%e.g [a,b]= dealarray([1,2]) assigns a=1, b=2

s = numel(arr);
n = nargout;

if n > s
    error('Insufficient number of elements in array!');
elseif n == 0
    return;
end

for i = 1:n
    varargout(i) = {arr(i)}; %#ok<AGROW>
end
end



