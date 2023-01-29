classdef TA_Plate_HX_m <TA_HX
    %Describes a parralel plate HX
    properties(AbortSet=true)
        %% independent properties
        Y0 (1,1) double{mustBeNonnegative}    % half the plate spacing  [m]
        Gas_Area_Ratio ...  % percentege of area covered by gas
            (1,1) double{mustBeGreaterThan(Gas_Area_Ratio,0)...
            ,mustBeLessThan(Gas_Area_Ratio,1.0001)}=0.1        
        %% dependent geometrical properties
        Diameter (1,1) double{mustBeNonnegative}  % HX diameter         [m]
        Area=0.1            % HX area [m^2]         [m]
        A_Solid=0.1         % Area covered by solid    [m]
        Lp (1,1) double{mustBeNonnegative}    % half thickness of solid [m]
    end
    
    methods
        %% constructor
        function obj = TA_Plate_HX_m(name,system,length,D,y0,phi,solid,Q)
            obj@TA_HX(name,system,length,solid,Q)
            obj.Diameter=D;
            obj.Y0=y0;
            obj.Gas_Area_Ratio=phi;
            obj.Lp=y0*(1/phi-1);
            obj.Area=pi*D^2/4;
            obj.A_Solid=pi*D^2/4*(1-obj.Gas_Area_Ratio);
        end
         %% set method for diameter
        function set.Diameter(obj,value)
            obj.Diameter=value;
            obj.Area=pi*value^2/4;
            obj.A_Solid=(1-obj.Gas_Area_Ratio)*obj.Area;
        end
        %% set method for area
        function set.Area(obj,value)
            obj.Area=value;
            obj.Diameter=(value/pi*4)^0.5;
            obj.A_Solid=(1-obj.Gas_Area_Ratio)*value;
        end
        %% set method for A_solid
        function  set.A_Solid(obj,value)
            if value>=obj.Area
                error('solid area must be smaller than total area')
            end
            obj.A_Solid=value;
            obj.Gas_Area_Ratio=(obj.Area-value)/obj.Area;
        end
        %% set method for Gas_Area ratio
        function  set.Gas_Area_Ratio(obj,value)
            if value>=1
                error('solid area must be smaller than total area')
            end
            obj.Gas_Area_Ratio=value;
            obj.A_Solid=(1-value)*obj.Area;
        end
        %% geometrical functions
        %this function  assists in averaging a property based
        %on the thermal penetration depth and geometry
        %source:delta E.C version 6.4b2.7 user's guide, equations
        %10.59-10.60

        function  f=f_function(obj,delta)
        rh=obj.Y0;
        f=tanh((1+1i)*rh/delta)/((1+1i)*rh/delta);
        end
        function epsilon=eps_function(obj,deltak,deltas,~)
            rh=obj.Y0;
            lp=obj.Lp;
            epsilon=tanh((1+1i)*rh/deltak)/tanh((1+1i)*lp/deltas);
         
        end
        %% solid temperature calculation
        function T_solid=Calculate_Solid_Temperature(obj)
          T_solid=0;
        end
    end
end

