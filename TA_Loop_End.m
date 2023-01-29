classdef  TA_Loop_End<TA_rget
    % marks a closing of a loop to go back to the beginning
    %allows targeting the difference between pressure, velocity and
    %temperature at the beginning and end.
    properties (SetAccess=protected)
        % the mass flux                    [W]
        Target_Value        %inverse of the acoustic impedance [-]
        
        Target_Definition=['The difference between pressure, velocity and'...
            'temperature at location and at begining of the system'...
            'usually targeted to be 0 for looped systems']
        IMRE    (1,1) logical=0
        Target_T=nan
        Target_P=nan
        Target_U=nan
    end
    
    methods
        %% constructor
        function obj = TA_Loop_End(name,system,IMRE)
            obj@TA_rget(name,system)
            obj.IMRE=IMRE;
        end
        %% calculation method
        function Calculate_Target(obj,~)
            if obj.IMRE
                obj.Target_Value=[obj.Pressure,obj.Velocity,obj.Temperature]-obj.System_H.Begin;
            else
                obj.Target_Value=[obj.Pressure,obj.Velocity,obj.Temperature]./...
                    obj.System_H.Begin;                
            end
            obj.Target_P=obj.Target_Value(1);
            obj.Target_U=obj.Target_Value(2);
            obj.Target_T=obj.Target_Value(3);
        end
    end
end

