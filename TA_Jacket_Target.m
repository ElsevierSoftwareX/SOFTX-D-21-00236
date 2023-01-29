classdef  TA_Jacket_Target<TA_rget
    % marks the beginning of a water jacket in the system.
    %allows to target the acoustic power and total power to be equal
    
    properties (SetAccess=protected)
        Target_Value
        Target_Definition=['Difference between total and acoustic power.'...
            'allows to target them to be equal when beginning a jacket']
    end
    
    methods
        %% constructor
        function obj = TA_Jacket_Target(name,system)
            obj@TA_rget(name,system)
        end
        %% calculation method
        function Calculate_Target(obj,H_before)
            if obj.Location==1
                error('This target cannot be the first component of a system')
            end
            P=obj.System_H.Components_H{obj.Location-1}.Pressure(end);
            U=obj.System_H.Components_H{obj.Location-1}.Velocity(end);
            E_before=0.5*real(P*conj(U));
            obj.Target_Value=H_before-E_before;
                
        end
        
    end
    
end

