classdef TA_Regenerator <TA_Component
    % Describes a dry regenerator in a TA system
    %valid only for a dry system. wet derivation is yet to be
    %performed
    properties
        %% independent properties 
        Length              % Stack length                      [m]
        Solid   (1,:)string % string containing solid name
        N_Sol_Points        % number of points in solution
    
        %% independent properties
        Area    (1,1) double ...                                      [m^2]
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}=0.1;                               % Regenerator area      [m^2]
        A_Solid(1,1) double{mustBeNonnegative}% Area covered by solid [m]    
        Rh   (1,1) double ...% duct hydraulic radius [m]
            {mustBeGreaterThan(Rh,1e-10),mustBeLessThan(Rh,1)}=0.0001;  
        K      (1,1) double{mustBeNonnegative}% porosity factor       [-]
    end
    
    methods
        %% constructor
        function obj = TA_Regenerator(name,system,length,A,A_solid,rh,K,solid)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Solid=solid;
            obj.Area=A;
            obj.A_Solid=A_solid;
            obj.Rh=rh;
            obj.K=K;
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            % begining as an input, and the variables at the end as an output
            
            %ode domain
            if isempty(obj.N_Sol_Points)
                domain=[locstart,locstart+obj.Length];
            else
                domain=locstart:obj.Length/(obj.N_Sol_Points-1):locstart+obj.Length;
            end
            [L,P_U_T] = ode113A(@(L,P_U_T) obj.ODEfun(P_U_T,H_before),domain,P_U_T_in);
            obj.X=L;
            obj.Pressure=P_U_T(:,1);
            obj.Velocity=P_U_T(:,2);
            obj.Temperature=P_U_T(:,3);
            PUT_After=P_U_T(end,:);
            H_after=H_before;
        end
        %% ODE function
        function output=ODEfun(obj,input,H_2)
            %this function integrates a Regenerator, including the effects of
            % viscocity and the thermoacoustic effect
            if ~obj.System_H.Dry_Switch
                error(['Regenerators can only',...
                     ' be added to dry systems at the moment'])
            end
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            Array=obj.System_H.Mixture_Array;
            
            %collect variables from regenerator
            A=obj.Area;
            A_solid=obj.A_Solid;
            solid =obj.Solid;
            rh=obj.Rh;
            ks_frac=obj.K;
            %calculate properties -mixture and solid
            [ nu, ~, ~, rho, Pr, ~, gamma, ~, ~,cp,k,sp] =...
                Mixture_Properties( p_m, abs(input(3)), 1,Array);
            mu=nu*rho;

            [k_solid,c_solid,rho_solid] =...
                solid_properties(solid,abs(input(3)));
            beta=1/abs(input(3));%thermal expansion coefficient =1/T
            
            
            u_1=input(2)/(A-A_solid);     %velocity
            phi=(A-A_solid)/A;            %porosity
            
            c1=1268-3545*phi+2544*phi^2;   %constants for calculation
            c2=-2.82+10.7*phi-8.6*phi^2;
            bphi=3.81-11.29*phi+9.47*phi^2;
            
            NR=4*abs(u_1)*rh*rho/mu;       %reynolds number
            epsilon_s=phi*rho*cp/(1-phi)/rho_solid/c_solid;
            epsilon_h=8*1i*rh^2/bphi/Pr^(1/3)/(2*k/omega/rho/cp);
            theta_P=angle(u_1)-angle(input(1));
            
            %integrals

             z=0:0.01:pi/2;
             g_c= 2/pi*trapz(z,1./(1+NR^0.6*cos(z).^0.6));    
             g_v= -2/pi*trapz(z,cos(2*z)./(1+NR^0.6*cos(z).^0.6));
            
            theta_T=angle(u_1);
            for i=1:5
            alpha1=(epsilon_s+(g_c+exp(2*1i*theta_P)*g_v)*epsilon_h)/(1+epsilon_s+epsilon_h*(g_c+exp(2*1i*theta_T)*g_v));
            alpha2=(epsilon_s+(g_c-g_v)*epsilon_h)/(1+epsilon_s+epsilon_h*(g_c+exp(2*1i*theta_T)*g_v));
            F(1)=-1i*omega*rho*u_1*(1+(1-phi)^2/2/(2*phi-1))-mu/rh^2*(c1/8+c2*NR/3/pi)*u_1;
            F(3)=(real((input(3)*beta*alpha1+1-input(3)*beta)*input(1)*conj(u_1))-2*H_2/(A-A_solid)) / (rho*cp/omega*imag(alpha2)*(u_1*conj(u_1))+2*(ks_frac*k_solid*(1-phi)/phi+k));
            TUM=(abs(input(3))*beta/rho/cp*alpha1*input(1)-alpha2/1i/omega*F(3)*u_1);
            F(2)=(-1i*omega*gamma/rho/sp^2*input(1)+beta*F(3)*u_1+1i*omega*beta*TUM)*(A-A_solid);%
            theta_T=angle(u_1)-angle(TUM);
            end
            
            
            output=F(:);
        end
        %% calculate derived properties
        function Calculate_Derived(obj)
            if isempty(obj.Pressure)
                error('system must be run before data is collected')
            end
            obj.Acoustic_Power=0.5*real(obj.Pressure.*conj(obj.Velocity));
            if obj.Location>1
                obj.Total_Power=obj.System_H.Components_H{obj.Location-1}.Total_Power(end)*ones(length(obj.Acoustic_Power),1);
            else
                obj.Total_Power=obj.Acoustic_Power(1)*ones(size(obj.Pressure));
            end
            % mass flux is 0
            obj.Mass_Flux=0*obj.Pressure;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a mesh stack, gas_area = Area - Solid_Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*(obj.Area-obj.A_Solid);
        end
    end
end

