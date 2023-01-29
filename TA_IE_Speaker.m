classdef  TA_IE_Speaker<TA_Component
    %abstract class, describes an thermoacoustic component
    %% constant properties
    properties
        Length=0            % length                                 [m]
        N_Sol_Points=2      % number of points in solution
    
        %% unprotected properties
        
        % area  [m^2]
        Area(1,1) double...
            {mustBeGreaterThan(Area,1e-10)}=0.1 
        
        % Coil electric resistence [ohm]
        R_E(1,1) double{mustBeNonnegative} 
        
        % Coil electric inertence  [h]
        L  (1,1) double{mustBeNonnegative} 
        
        % maximal force  unit current  [N/A]
        % this is the force produced by the 
        % coil when it is held firmly in place.
        F_over_I (1,1) double{mustBeNonnegative}            
        
        % moving mass of the transducer [kg]
        M (1,1) double{mustBeNonnegative}    
        
        % spring constant of the transducer      [N/m]
        K (1,1) double{mustBeNonnegative}
        
        % mechanical resistence of  transducer   [Ns/m]
        R_M (1,1) double{mustBeNonnegative}    
        
        % complex current input on loudspeaker   [A]
        I   (1,1) double    
        V   (1,1) double =nan;
        %electric consumption                    [W]
        Electric_Consumption (1,1) double
        
        Jacket  (1,1) logical   %whether or not the loudspeaker is Jacketed
        
    end
    properties (Hidden)
        % an assisting property for using insertsmooth and removesmooth
        %should be removed once links are available
        %when equal to 1 the loudspeaker is "silenced" and has no effect.
        %when equal to 0 the loudspeaker functions regularly.
        Silencer (1,1) double=0;
    end
    methods
        function obj=TA_IE_Speaker(name,system,Area,R_E,L,F_over_I,M,k,R_M,I)
            obj@TA_Component(name,system)
            obj.Area=Area;
            obj.R_E=R_E;
            obj.L=L;
            obj.F_over_I=F_over_I;
            obj.M=M;
            obj.K=k;
            obj.R_M=R_M;
            obj.I=I;      
        end
          %% set method- to keep length and nsolpoint constant
        function set.N_Sol_Points(~,~)
            error('N_Sol_Points for a Speaker must be 2')
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            %source-delta E.C version 6.4b2.7 user's guide, section 10.3.1
            %collect variables from system
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            dryswitch=obj.System_H.Dry_Switch;
            Array=obj.System_H.Mixture_Array;
            
            %collect variables from object
            A=obj.Area;
            T=P_U_T_in(3);
            Bl=obj.F_over_I;
            I1=obj.I;
            Re=obj.R_E;
            Rm=obj.R_M;
            
            %mixture properties
            [ ~, alpha, ~, rho, ~, ~, gamma, ~, ~,~,~,sp] = Mixture_Properties(p_m,T,dryswitch,Array);
            
            deltak=(2*alpha/omega)^0.5;
            tau=Bl/A;
            tau_prime=-tau;
            Ze=Re+1i*omega*obj.L;
            Zm=Rm/A^2+1i*(omega*obj.M-obj.K/omega)/A^2;
            
            U1=P_U_T_in(2)-omega/rho/sp^2*(gamma-1)*deltak/2*A*P_U_T_in(1);
            V1=Ze*I1-tau*U1;
            obj.V=V1;
            p1_out=P_U_T_in(1)+(tau_prime*I1-Zm*U1)*(1-obj.Silencer);
            U1_out=P_U_T_in(2)-omega/rho/sp^2*(gamma-1)*deltak/2*A*(P_U_T_in(1)+p1_out)*(1-obj.Silencer);
            obj.Electric_Consumption=0.5*real(I1*conj(V1));
            
            
            obj.X=[locstart,locstart+obj.Length];
            obj.Pressure=[P_U_T_in(1);p1_out*ones(1,obj.N_Sol_Points-1)];
            obj.Velocity=[P_U_T_in(2);U1_out*ones(1,obj.N_Sol_Points-1)];
            obj.Temperature=P_U_T_in(3)*ones(obj.N_Sol_Points,1);
            PUT_After=[p1_out,U1_out,T];
            if obj.Jacket
            H_after=0.5*real(p1_out*conj(U1_out));
            else
                H_after=H_before+obj.Electric_Consumption;
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
            else 
                if obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
                else
                TP=obj.Acoustic_Power(1);
                end
                obj.Total_Power=[TP;TP+obj.Electric_Consumption];
            end
            obj.Mass_Flux=0*obj.X;   
        end
        %% addition to Empty variables for Speakers
        function Empty_Varibles(obj)
            obj.Empty_Varibles@TA_Component;
            obj.Electric_Consumption=0;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a loudspeaker, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
    end
    
end

