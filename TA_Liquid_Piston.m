classdef  TA_Liquid_Piston<TA_Component
    %compliance   
    properties
        %% unprotected properties
        Length=0            % component length- only for plotting   [m]
        N_Sol_Points=2      % number of points in solution
        Rho   (1,1) double{mustBeNonnegative} % liquid density   [kg/m^3]
        Mu    (1,1) double{mustBeNonnegative} % dynamic viscosity [pa s]
        Area  (1,1) double{mustBeNonnegative} % Piston area       [m^2]
        Ltotal ...                           %total water collumn length[m]  
           (1,1) double{mustBeNonnegative}        
        Jacket(1,1) logical  % whether or not the liquid piston is jacketed       
    end
       
    methods
        %% constructor
        function obj=TA_Liquid_Piston(name,system,length,rho,mu,A,L)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Rho=rho;
            obj.Mu=mu;
            obj.Area=A;
            obj.Ltotal=L;
        end
        %% set method to keep nsolpoints constant
          function set.N_Sol_Points(~,~)
            error('N_Sol_Points for a liquid piston must be 2')
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            
            %collect variables from system
            omega=obj.System_H.Frequency*2*pi;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %collect variables from object
            rho=obj.Rho;
            mu=obj.Mu;
            A=obj.Area;
            L=obj.Ltotal;
            
            
            D=(A*(4/pi))^0.5;
            E=0.25;%E values associated with the total minor losses in the U-shaped tube are estimated to be 0.25.
            r_t=D/2;%tube radius [m]
            Freq=omega/(2*pi);
            X_l=abs(P_U_T_in(2))/(omega*A);%liquid displacemnt
            V_l=A*L;%liquid volume
            g=9.8;%gravity[m/s^2]
            
            R_v=(pi^1.5)*D*L*(rho*mu*Freq)^0.5;
            R_k=0.84*pi^2*rho*Freq*X_l*r_t^2*E;
            R_m=R_v+R_k;
            
            P_after=1*P_U_T_in(1)+(-1/A^2)*(R_m+1i*(omega*rho*V_l-(2*rho*A*g)/omega))*P_U_T_in(2); 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            PUT_After=[P_after,P_U_T_in(2),P_U_T_in(3)];
            obj.X=[locstart,locstart+obj.Length];
            obj.Pressure=[P_U_T_in(1);P_after*ones(1,obj.N_Sol_Points-1)];
            obj.Velocity=P_U_T_in(2)*ones(obj.N_Sol_Points,1);
            obj.Temperature=P_U_T_in(3)*ones(obj.N_Sol_Points,1);
            
            
            if obj.Jacket
                H_after=0.5*real(conj(PUT_After(2))...
                    *PUT_After(1));
            else
                H_after=H_before;
            end
        end
        %% calculate derived properties
        %this function calculates the derived properties (mass flux,
        %acoustic power, total power) based on the pressure, temperature
        % and additional parameters
        function Calculate_Derived(obj)
            if isempty(obj.Pressure)
                error('system must be run before data is collected')
            end
            obj.Acoustic_Power=0.5*real(obj.Pressure.*conj(obj.Velocity));
            
            %calculate total power
            if obj.Jacket
                obj.Total_Power=obj.Acoustic_Power;
                if obj.Location>1&&obj.Total_Power(1) ~=obj.System_H.Components_H{obj.Location-1}.Total_Power(end)
                    warning('jacket induces jump in acoustic power. consider removing or adding a jacktarget')        
                end
            elseif obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
                obj.Total_Power=TP*ones(size(obj.Pressure));
            else
                obj.Total_Power=obj.Acoustic_Power(1)*ones(size(obj.Pressure));
            end
            obj.Mass_Flux=0*obj.X;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a liquid piston, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
    end
    
end

