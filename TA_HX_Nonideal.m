classdef TA_HX_Nonideal <TA_Component
    %Abstract component, describes a heat exchanger in a TA system
    % subclasses:TA_Plate_HX
    properties (AbortSet=true)
        %% independent properties
        Length                         % Stack length                   [m]
        Heat_Input (1,1)double         % Heat input into heat exchanger [W]
        SolidT...                      % temperature of the solid[k]
            (1,1) double{mustBeNonnegative}
        Solid      (1,:)string         % string containing solid name
        N_Sol_Points                   % number of points in solution
        
    end
    properties(Abstract)
        %% Abstract properties
        Area (1,1) double ...                % HX area            [m^2]
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}  
        
        A_Solid   (1,1) double ...           %HX solid Area       [m]
            {mustBeNonnegative,mustBeLessThan(A_Solid,1)}  
    end
    
    properties (Hidden)
        %% Hidden properties for solid_properties_ODE
        saved_s_properties=NaN(3,1);
        saved_T_s = NaN;
    end
    
    methods
        %% constructor
        function obj = TA_HX_Nonideal(name,system,length,solid,heat)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Solid=solid;
            obj.Heat_Input=heat;
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            
            %ode domain
            if isempty(obj.N_Sol_Points)
                domain=[locstart,locstart+obj.Length];
            else
                domain=locstart:obj.Length/(obj.N_Sol_Points-1):locstart+obj.Length;
            end
            [L,P_U_T] = ode113A(@(L,P_U_T) obj.ODEfun(P_U_T),domain,P_U_T_in);
            obj.X=L;
            obj.Pressure=P_U_T(:,1);
            obj.Velocity=P_U_T(:,2);
            obj.Temperature=P_U_T(:,3);
            obj.SolidT=obj.Calculate_Solid_Temperature;
            PUT_After=P_U_T(end,:);
            H_after=H_before+obj.Heat_Input;
        end
        %% ODE function
        function output=ODEfun(obj,input)
            %this function integrates a Heat exchanger, including the effects of
            % viscocity and the thermoacoustic effect
            
            %source
            %   Theoretical performance characteristics of a travelling-
            %      wave phase-change thermoacoustic engine for low-grade
            %      heat recovery, Rui Yang,Avishai Meir, and Guy Ramon
            %      Applied energy 2020
            
            
            %input parameter- a vetor containing the complex pressure and
            %velocity and the mean temperature at a point
            %input(1)--p1,input(2)--U1,input(3)--T_m
            %output paremeter-the values of dp/dx, dT/dx and dU/dx at this
            %point
            % collect variables from system
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            gas=obj.System_H.Inert;

            %collect variables from stack
            A=obj.Area;
            A_solid=obj.A_Solid;
            solid =obj.Solid;
            
            %calculate properties -mixture and solid
             rho=Modrefpropm('D','T',abs(input(3)),'P',p_m/1000,gas); %kg/m3
            cp=Modrefpropm('C','T',abs(input(3)),'P',p_m/1000,gas);%J/kg/k
            nu=0.0001*Modrefpropm('$','T',abs(input(3)),'P',p_m/1000,gas); %m2/s
            k=Modrefpropm('L','T',abs(input(3)),'P',p_m/1000,gas);
            alpha=0.0001*Modrefpropm('%','T',abs(input(3)),'P',p_m/1000,gas);%m2/s
            sp=Modrefpropm('A','T',abs(input(3)),'P',p_m/1000,gas);
            gamma=Modrefpropm('K','T',abs(input(3)),'P',p_m/1000,gas);

            [k_solid,c_solid,rho_solid] = obj.solid_properties_ODE(solid,abs(input(3)));
            
            deltanu=(2*nu/omega)^0.5;
            deltak=(2*alpha/omega)^0.5;

            fnu=obj.f_function(deltanu);
            fk=obj.f_function(deltak);
            
            %fnu=tanh((1+1i)*rh/deltanu)/((1+1i)*rh/deltanu);
            %fk=tanh((1+1i)*rh/deltak)/((1+1i)*rh/deltak);
            
            
            delta_s=(2*k_solid/omega/rho_solid/c_solid)^0.5;
            epsilon_s=(k*rho*cp/k_solid/rho_solid/c_solid)^0.5*eps_function(obj,deltak,delta_s,fk);
            

            %F(1)--dP1/dx,F(2)--dU1/dx,F(3)=dTm/dx
            %equations 1-3 without a temperature gradient
            F(1)=-1i*omega*rho*input(2)/(1-fnu)/(A-A_solid);
            F(2)=-1i*omega*(A-A_solid)/rho/sp^2*(1+(gamma-1)*fk/(1+epsilon_s))*input(1);
            F(3)=0;
            output=F(:);
        end
        %% Solid properties for ODE
        function [k,cp,rho] = solid_properties_ODE(obj,solid,T)
            if T == obj.saved_T_s
                k=obj.saved_s_properties(1); cp=obj.saved_s_properties(2);
                rho=obj.saved_s_properties(3);
            else
                [k,cp,rho] = solid_properties(solid,T);
                obj.saved_T_s=T;
                obj.saved_s_properties=[k,cp,rho];
            end
        end

        % initializing saved properties when Solid is modified
        function set.Solid(obj,value)
            obj.Solid=value;
            obj.saved_T_s=NaN;
            obj.saved_s_properties=NaN(3,1);
        end
        %% calculate derived properties
        function Calculate_Derived(obj)
            if isempty(obj.Pressure)
                error('system must be run before data is collected')
            end
            obj.Acoustic_Power=0.5*real(obj.Pressure.*conj(obj.Velocity));
            
            %total power changes linearly between beginning and end
            if obj.Location>1
                Tp1=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
            else
                Tp1=obj.Acoustic_Power(1);
            end
            Tp2=Tp1+obj.Heat_Input;
            obj.Total_Power=linspace(Tp1,Tp2,length(obj.Acoustic_Power))';
            %caculate mass flux
            obj.Mass_Flux=0*obj.Pressure;
            
            %internal function
            
        end
        %% addition to Empty variables for HXs
        function Empty_Varibles(obj)
            obj.Empty_Varibles@TA_Component;
            obj.SolidT=0;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a HX, gas_area = Area - Solid_Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*(obj.Area-obj.A_Solid);
        end
    end
    methods(Abstract)
        %each stack must have functions to represent the geometry of the
        % stack, in term of the parameter f and epsilon.
        % see-delta E.C version 6.4b2.7 user's guide, section 10.5.1
        
        %f function, representing spactial average inside pore
        f=f_function(obj,delta)
        %epsilon function, representing solid activity divided by the
        %(always identical) first term.
        e=eps_function(obj,deltak,deltas,fk)
        
        T=Calculate_Solid_Temperature(obj)
    end
    
end

