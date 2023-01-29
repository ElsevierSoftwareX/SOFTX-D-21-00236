classdef  TA_Loop_End_PHAB<TA_rget
    % marks a closing of a loop to go back to the beginning
    %allows targeting the difference between pressure, velocity and
    %temperature at the beginning and end.
    properties (SetAccess=protected)
        % the mass flux                    [W]
        Target_Value        %inverse of the acoustic impedance [-]
        
        Target_Definition=['The difference between pressure, velocity and'...
            'temperature at location and at begining of the system'...
            'usually targeted to be 1 for looped systems']                             
    end
    
    methods
        %% constructor
        function obj = TA_Loop_End_PHAB(name,system)
            obj@TA_rget(name,system)
        end
        %% calculation method
        function Calculate_Target(obj,~)
           obj.Target_Value=[obj.Pressure,obj.Velocity,obj.Temperature]./...
               obj.System_H.Begin;
        end     
    end     
end

