classdef TA_Yuval_Practice < TA_Component
    % Practice component which add 1 to all values

    properties
        %%
        Length=0                 % Duct length [m]
        N_Sol_Points=0           % number of points in solution
        Mood {mustBeA(Mood,["cell","string","char"])} = 'alright' % The mood of the component
    end

    methods
        %% constructor
        function obj = TA_Yuval_Practice(name,system,M)
            obj@TA_Component(name,system)
            obj.Mood = M;
        end

        %% main method
        function [PUT_After,H_after] = run_component(obj,P_U_T_in,H_in,locstart)
            obj.Length = obj.Length + 1;
            obj.N_Sol_Points = obj.N_Sol_Points + 1;
            obj.X = obj.X + 1;
            obj.Pressure = obj.Pressure + 1;
            obj.Velocity = obj.Velocity + 1;
            obj.Temperature = obj.Temperature + 1;
            obj.Acoustic_Power = obj.Acoustic_Power +1;
            obj.Mass_Flux = obj.Mass_Flux + 1;
            obj.Total_Power = obj.Total_Power + 1;
            obj.Mood = "very happy";
            % adding outputs for the runSystem to work
            PUT_After = P_U_T_in(end,:);
            H_after = H_in(end,:);
        end

        % remaining abstract methods for the component to work
        
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