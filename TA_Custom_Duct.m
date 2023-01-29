classdef TA_Custom_Duct <TA_Duct
    %describes a duct in a thermoacoustic system
    properties (AbortSet=true)
        %% dependent geometrical properties
        Area=0.1            % duct area                 [m^2]
        Rh=0.1              % duct hydraulic radius     [m]
    end   
    
    methods
        %% constructor
        function obj = TA_Custom_Duct(name,system,length,area,rh)
            obj@TA_Duct(name,system,length)
            obj.Area=area;
            obj.Rh=rh;
        end
    end  
end

