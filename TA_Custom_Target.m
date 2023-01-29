classdef  TA_Custom_Target<TA_rget
    % marks the beginning of a water jacket in the system.
    %allows to target the acoustic power and total power to be equal
    
    properties (SetAccess=protected)
        Target_Value 
        Target_Definition='Custom Target defined by the inline function f'
        F (1,1) function_handle=@()1+1
    end
    
    methods
        %% constructor
        function obj = TA_Custom_Target(name,system,f)
            obj@TA_rget(name,system)
            obj.F=f;
        end
        %% calculation method
        function Calculate_Target(obj,~)
            obj.Target_Value=obj.F();
        end
        
    end
    
end

