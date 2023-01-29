classdef  TA_Compliance<TA_Component
    %compliance   
    properties
        %% unprotected properties
        Length=0            % Compliance length- only for plotting   [m]
        N_Sol_Points=2      % number of points in solution
        Volume  (1,1) double{mustBeNonnegative}% compliance volume   [m^3]
        
        Area  (1,1) double{mustBeNonnegative} % compliance area       [m^2]
                   
        Jacket(1,1) logical     % whether or not the compliance is jacketed       
    end
       
    methods
        %% constructor
        function obj=TA_Compliance(name,system,V,A)
             obj@TA_Component(name,system)
            obj.Area=A;
            obj.Volume=V;
        end
        %% set method to keep nsolpoints constant
          function set.N_Sol_Points(~,val)
              if val ~= 2
                error('N_Sol_Points for a compliacne must be 2')
              end
          end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            
            %collect variables from system
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            dryswitch=obj.System_H.Dry_Switch;
            Array=obj.System_H.Mixture_Array;
            
            %collect variables from object
            S_cp=obj.Area;
            V_cp=obj.Volume;
            T=P_U_T_in(3);
            
            %calculation- source:
            % GST- Gregory swift-Thermoacoustics, A Unifying Perspective,
            %      second edition, chapter 4.1 Springer 2017
            
            %source-delta E.C version 6.4b2.7 user's guide, section 10.2.1
            [ ~, alpha, ~, rho, ~, ~, gamma, ~, ~,~,~,sp] = Mixture_Properties(p_m,T,dryswitch,Array);
            deltak=(2*alpha/omega)^0.5;
            
            U1_after=P_U_T_in(2)-1i*omega*P_U_T_in(1)/rho/sp^2*(V_cp-1i*(gamma-1)*S_cp*deltak/2);
            PUT_After=[P_U_T_in(1),U1_after,P_U_T_in(3)];
            obj.X=[locstart,locstart+obj.Length];
            obj.Pressure=P_U_T_in(1)*ones(obj.N_Sol_Points,1);
            obj.Velocity=[P_U_T_in(2);U1_after*ones(1,obj.N_Sol_Points-1)];
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
        %In this case of a compliance, we make the assumption that gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
    end
    
end

