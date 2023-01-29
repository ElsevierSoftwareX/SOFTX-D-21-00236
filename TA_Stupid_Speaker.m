classdef  TA_Stupid_Speaker<TA_Component
    %abstract class, describes an thermoacoustic component
    %% constant properties
    properties
        Length=0            % length                                 [m]
        N_Sol_Points=2      % number of points in solution
    
        %% unprotected properties
        
        % area  [m^2]
        Area(1,1) double...
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}=0.1 
        
        % Coil electric resistence [ohm]
        Z_m(1,1) double
        I
        Tau
        Jacket  (1,1) logical   %whether or not the loudspeaker is Jacketed
        V
        L
        R
    end
    
    methods
        function obj=TA_Stupid_Speaker(name,system,Z_m,I,Tau,L,R)
            obj@TA_Component(name,system)
            obj.Z_m=Z_m;
            obj.I=I;
            obj.Tau=Tau;
            obj.L=L;
            obj.R=R;
        end
          %% set method- to keep length and nsolpoint constant
        function set.Length(~,~)
            error('Length of a Speaker must be zero')
        end
         function set.N_Sol_Points(~,~)
            error('N_Sol_Points for a Speaker must be 2')
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            %source-delta E.C version 6.4b2.7 user's guide, section 10.3.1
            %collect variables from system

            %collect variables from object
            T=P_U_T_in(3);

            
            %mixture properties            
            tau=obj.Tau;
            tau_prime=-tau;
            Zm=obj.Z_m;
            
            U1=P_U_T_in(2);
            I1=obj.I;
            p1_out=P_U_T_in(1)+tau_prime*I1-Zm*U1;
            U1_out=P_U_T_in(2);
            
            
            obj.X=locstart*ones(obj.N_Sol_Points,1);
            obj.Pressure=[P_U_T_in(1);p1_out*ones(1,obj.N_Sol_Points-1)];
            obj.Velocity=[P_U_T_in(2);U1_out*ones(1,obj.N_Sol_Points-1)];
            obj.Temperature=P_U_T_in(3)*ones(obj.N_Sol_Points,1);
            PUT_After=[p1_out,U1_out,T];
            if obj.Jacket
            H_after=0.5*real(p1_out*conj(U1_out));
            else
                H_after=H_before;
            end
            omega=obj.System_H.Frequency*2*pi;
            Ze=obj.R+1i*omega*obj.L;
            obj.V=Ze*obj.I-tau*U1;
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
                obj.Total_Power=[TP;TP];
            end
            obj.Mass_Flux=0*obj.X;   
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a speaker, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
    end
    
end

