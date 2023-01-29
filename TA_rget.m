classdef  TA_rget<TA_Component
    % any type of target in a TA system
    
    properties
        %% unprotected properties
        Length=0         % length [m]
        N_Sol_Points=1   % number of points in solution
        
    end
    properties (Abstract,SetAccess=protected)
        Target_Value    double     %inverse of the acoustic impedance [-]
        Target_Definition (1,:) char   % string containing target definition
    end
    
    methods
        function obj = TA_rget(name,system)
            obj@TA_Component(name,system)
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            obj.X=locstart*ones(length(obj.N_Sol_Points),1);
            obj.Pressure=P_U_T_in(1);
            obj.Velocity=P_U_T_in(2);
            obj.Temperature=P_U_T_in(3);
            PUT_After=P_U_T_in;
            H_after=H_before;
            try
            obj.Calculate_Target(H_before);
            catch me
            obj.Pressure=[];
            obj.Velocity=[];
            obj.Temperature=[];
            rethrow(me)
            end
        end
        %%
        function set.Length(~,val)
            if val ~= 0
                error('Length of a target must be zero')
            end
        end
        function set.N_Sol_Points(~,val)
            if val ~= 1
                error('N_Sol_Points for a target must be one')
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
            if obj.Location>1
                obj.Total_Power=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
            else
                obj.Total_Power=obj.Acoustic_Power(1);
            end
            obj.Mass_Flux=0*obj.X;
        end
        %% addition to Empty variables for Targets
        function Empty_Varibles(obj)
            obj.Empty_Varibles@TA_Component;
            obj.Target_Value=0;
        end
        %% get gas area
        function gas_area=get_gas_area(obj)
            gas_area=NaN;
        end
    end
    methods (Abstract)
        Calculate_Target(obj,H)
    end
end

