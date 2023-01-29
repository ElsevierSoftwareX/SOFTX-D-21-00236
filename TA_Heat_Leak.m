classdef  TA_Heat_Leak<TA_Component
    %compliance   
    properties
        %% unprotected properties
        Length=0            % component length- only for plotting   [m]
        N_Sol_Points=2      % number of points in solution
        Area  (1,1) double{mustBeNonnegative} % Piston area       [m^2]
        H     (1,1) double{mustBeNonnegative} % heat transfer coefficient
        Q_Loss       
        T_Out (1,1)
    end
       
    methods
        %% constructor
        function obj=TA_Heat_Leak(name,system,A,h,T,L)
            obj@TA_Component(name,system)
            obj.Length=L;
            obj.Area=A;
            obj.H=h;
            obj.T_Out=T;
        end
        %% set method to keep nsolpoints constant
          function set.N_Sol_Points(~,~)
            error('N_Sol_Points for a liquid piston must be 2')
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output

            
            PUT_After=P_U_T_in;
            obj.X=[locstart,locstart+obj.Length];
            obj.Pressure=P_U_T_in(1)*ones(obj.N_Sol_Points,1);
            obj.Velocity=P_U_T_in(2)*ones(obj.N_Sol_Points,1);
            obj.Temperature=P_U_T_in(3)*ones(obj.N_Sol_Points,1);
            obj.Q_Loss=obj.H*obj.Area*(P_U_T_in(3)-obj.T_Out);
            H_after=H_before-obj.Q_Loss;
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
           if obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end)-obj.Q_Loss;
            else
                TP=obj.Acoustic_Power(1)-obj.Q_Loss;
           end
           obj.Total_Power=TP*ones(length(obj.Pressure),1);
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

