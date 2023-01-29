classdef TA_Linear_Cone <TA_Cone
    %describes a duct in a thermoacoustic system
    properties (AbortSet=true)
        %% dependent geometrical properties
        Diameter_1  (1,1) double{mustBeNonnegative} % first diameter  [m]  
        Diameter_2  (1,1) double{mustBeNonnegative} % second diameter [m]
        Area_1      (1,1) double{mustBeNonnegative} % first area     [m^2]
        Area_2      (1,1) double{mustBeNonnegative} % second area    [m^2]
        Rh_1        (1,1) double{mustBeNonnegative}  % hydraulic radiuses [m]
        Rh_2        (1,1) double{mustBeNonnegative}        
    end   
    
    methods
        %% constructor
        function obj = TA_Linear_Cone(name,system,length,D1,D2)
            obj@TA_Cone(name,system,length)
            obj.Diameter_1=D1;
            obj.Diameter_2=D2;
        end
        %% set method for diameter
        function set.Diameter_1(obj,value)
            obj.Diameter_1=value;
            obj.Rh_1=value/4;
            obj.Area_1=pi*value^2/4;
        end
          %% set method for diameter
        function set.Diameter_2(obj,value)
            obj.Diameter_2=value;
            obj.Rh_2=value/4;
            obj.Area_2=pi*value^2/4;
        end
        %% set method for area
        function set.Area_1(obj,value)
            obj.Area_1=value;
            obj.Rh_1=(value/pi)^0.5/2;
            obj.Diameter_1=(value/pi*4)^0.5;
        end
        %% set method for area
        function set.Area_2(obj,value)
            obj.Area_2=value;
            obj.Rh_2=(value/pi)^0.5/2;
            obj.Diameter_2=(value/pi*4)^0.5;
        end
        %% set method for rh
        function  set.Rh_1(obj,value)
            obj.Rh_1=value;
            obj.Diameter_1=value*4;
            obj.Area_1=pi*value^2*4;         
        end
             %% set method for rh
        function  set.Rh_2(obj,value)
            obj.Rh_2=value;
            obj.Diameter_2=value*4;
            obj.Area_2=pi*value^2*4;         
        end
        %% calculation methods
        function A=Calc_A(obj,x_rel)
            r=obj.Rh_1+x_rel/obj.Length*(obj.Rh_2-obj.Rh_1);
            A=pi*r^2*4;
        end
        function rh=Calc_Rh(obj,x_rel)
            rh=obj.Rh_1+x_rel/obj.Length*(obj.Rh_2-obj.Rh_1);
            
        end
    end
end

