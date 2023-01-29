classdef  TA_Soft_End<TA_rget
    % marks a zero pressure at the end of the system
    %allows targeting the complex impedance at the wall to be 0
    %often more stable than targeting the pressure to be zero
    %more directly.
    %inverse of the acoustic impedance [-]
    properties (SetAccess=protected)
        Target_Value %(1,1)double
        Target_Definition=['The inverse of the impedance at location.'...
            'Usually used with 0 as target to represent a solid wall']
    end
    
    methods
        %% constructor
        function obj = TA_Soft_End(name,system)
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
                    'the Softend. Softend must follow some component with area'])
            end
            [ ~, ~, ~, rho, ~, ~, ~, ~, ~,~,~,sp] = Mixture_Properties(obj.System_H.P_m,obj.Temperature(end),obj.System_H.Dry_Switch,obj.System_H.Mixture_Array);          
            obj.Target_Value=1/obj.Velocity.*(obj.Pressure*Area)/rho/sp; %inverse of the acoustic impedance [-]
        end
        
    end
    
end

