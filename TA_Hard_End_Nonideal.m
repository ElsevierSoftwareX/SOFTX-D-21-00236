classdef  TA_Hard_End_Nonideal  <TA_rget
    % marks a solid wall at the end of the system
    %allows targeting the complex impedance at the wall to be 0
    %often more stable than targeting the velocity to be zero
    %more directly.
    %inverse of the acoustic impedance [-]
    properties (SetAccess=protected)
        Target_Value %(1,1)double
        Target_Definition=['The inverse of the impedance at location.'...
            'Usually used with 0 as target to represent a solid wall']
    end

    methods
        %% constructor
        function obj = TA_Hard_End_Nonideal(name,system)
            obj@TA_rget(name,system)
        end

        function Calculate_Target(obj,~)
            % find Area of last component with an area
            Area=0;
            for i=obj.Location-1:-1:1
                try
                    Area=obj.System_H.Components_H{i}.Area;
                    break
                catch
                    try
                        Area=obj.System_H.Components_H{i}.Area_2;
                        break
                    catch
                    end
                end
            end
            if Area==0
                error(['No component with an area has been found before',...
                    'the Hardend. Hardend must follow some component with area'])
            end
            gas=obj.System_H.Inert;
            rho=Modrefpropm('D','T',obj.Temperature(end),'P',obj.System_H.P_m/1000,gas); %kg/m3
            sp=Modrefpropm('A','T',obj.Temperature(end),'P',obj.System_H.P_m/1000,gas);
            obj.Target_Value=obj.Velocity./(obj.Pressure*Area)*rho*sp; %inverse of the acoustic impedance [-]
        end

    end

end

