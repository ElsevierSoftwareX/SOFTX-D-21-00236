classdef TA_Stack_Nonideal <TA_Component
    %abstract component, describes a stack in a TA system
    % subclasses:TA_Plate_Stack, TA_Rectangular_stack
    properties
        %% independent properties
        Length              % Stack length            [m]
        Solid  (1,:)string  % string containing solid name
        N_Sol_Points        % number of points in solution
    end
       properties(Abstract)
        %% Abstract properties
        Area (1,1) double ...                % Stack area            [m^2]
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}  
        
        A_Solid   (1,1) double ...           %Stack solid Area       [m]
            {mustBeNonnegative,mustBeLessThan(A_Solid,1)}  
    end
    
    methods
        %% constructor
        function obj = TA_Stack_Nonideal(name,system,length,solid)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Solid=solid;
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
            %this function integrates a Stack, including the effects of
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
            k=Modrefpropm('L','T',abs(input(3)),'P',p_m/1000,gas);
            cp=Modrefpropm('C','T',abs(input(3)),'P',p_m/1000,gas);%J/kg/k
            nu=0.0001*Modrefpropm('$','T',abs(input(3)),'P',p_m/1000,gas); %m2/s
            alpha=0.0001*Modrefpropm('%','T',abs(input(3)),'P',p_m/1000,gas);%m2/s
            sp=Modrefpropm('A','T',abs(input(3)),'P',p_m/1000,gas);
            gamma=Modrefpropm('K','T',abs(input(3)),'P',p_m/1000,gas);
            beta=Modrefpropm('B','T',abs(input(3)),'P',p_m/1000,gas);

            [k_solid,c_solid,rho_solid] = solid_properties(solid,abs(input(3)));
            
            deltanu=(2*nu/omega)^0.5;
            deltak=(2*alpha/omega)^0.5;
            Pr=(deltanu/deltak)^2;
            fnu=obj.f_function(deltanu);
            fk=obj.f_function(deltak);
            
            delta_s=(2*k_solid/omega/rho_solid/c_solid)^0.5;
            epsilon_s=(k*rho*cp/k_solid/rho_solid/c_solid)^0.5*eps_function(obj,deltak,delta_s,fk);

            %dry terms
            dry1=1+(gamma-1)*fk/(1+epsilon_s);
            dry2=beta*(fk-fnu)/(1-fnu)/(1-Pr)/(1+epsilon_s);
            dry3=0.5*real(input(1)*conj(input(2))*(1-(fk-conj(fnu))*beta*abs(input(3))/(1+epsilon_s)/(1+Pr)/(1-conj(fnu))));
            dry4=imag(conj(fnu)+(fk-conj(fnu))*(1+epsilon_s*fnu/fk)/(1+epsilon_s)/(1+Pr));

            %F(1)--dP1/dx,F(2)--dU1/dx,F(3)=dTm/dx
            %equations 1-3
            F(1)=-1i*omega*rho*input(2)/(1-fnu)/(A-A_solid);
           % F(3)=(-H_2+dry3) / (p_m/Rg/abs(input(3))*((abs(input(2)))^2/2/(A-A_solid)/omega/(abs(1-fnu))^2*(dry4+wet4)+wet5*(A-A_solid))+((A-A_solid)*k+A_solid*k_solid));
            F(3)=(H_2-dry3) / (rho*cp*abs(input(2))^2/2/omega/(A-A_solid)/(1-Pr)/abs(1-fnu)^2*(dry4)-((A-A_solid)*k+1*A_solid*k_solid));

            F(2)=-1i*omega*(A-A_solid)/rho/sp^2*(dry1)*input(1)+(dry2)*F(3)*input(2);%
            
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
            %caculate mass flux
            obj.Mass_Flux=0*obj.Pressure;
        end
        
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a stack, gas_area = Area - Solid_Area.
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
        %epsilon function, representing solid activity divided by the fluid
        %(always identical) first term.
        e=eps_function(obj,deltak,deltas,fk)
    end
    
end

