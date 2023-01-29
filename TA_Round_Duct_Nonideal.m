classdef TA_Round_Duct_Nonideal <TA_Duct_Nonideal
    %describes a duct in a thermoacoustic system
    properties (AbortSet=true,SetObservable=true)
        %% dependent geometrical properties
        Diameter...                                % duct diameter    [m]
            (1,1) double{mustBeNonnegative}=0.10032  
        Area=0.10032                               % duct area        [m^2]
        Rh=0.10032                                 % hydraulic radius [m]
    end   
    
    methods
        %% constructor
        function obj = TA_Round_Duct_Nonideal(name,system,length,diameter)
            obj@TA_Duct_Nonideal(name,system,length)
            obj.Diameter=diameter;
        end
        %% set method for diameter
        function set.Diameter(obj,value)
            obj.Diameter=value;
            obj.Rh=value/4;
            obj.Area=pi*value^2/4;
        end
        %% set method for area
        function set.Area(obj,value)
            obj.Area=value;
            obj.Rh=(value/pi)^0.5/2;
            obj.Diameter=(value/pi*4)^0.5;
            notify(obj,'Area_Changed');
        end
        %% set method for rh
        function  set.Rh(obj,value)
            obj.Rh=value;
            obj.Diameter=value*4;
            obj.Area=pi*value^2*4;         
        end
    end
    events
        Length_Changed
        Area_Changed
    end
end

