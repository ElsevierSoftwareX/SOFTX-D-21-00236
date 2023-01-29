classdef TA_Custom_Cone <TA_Cone
    %describes a duct in a thermoacoustic system
    properties (AbortSet=true)
        %% function handles
        Area_Function  (1,1) function_handle=@()1+1
        RH_Function  (1,1) function_handle=@()1+1

    end   
    
    methods
        %% constructor
        function obj = TA_Custom_Cone(name,system,length,fa,fr)
            obj@TA_Cone(name,system,length)
            obj.Area_Function=fa;
            obj.RH_Function=fr;
        end
        %% calculation methods
        function A=Calc_A(obj,x_rel)
            A=obj.Area_Function(x_rel);
        end
        function rh=Calc_Rh(obj,x_rel)
            rh=obj.RH_Function(x_rel);
            
        end
    end
end

