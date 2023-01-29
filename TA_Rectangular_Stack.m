classdef TA_Rectangular_Stack <TA_Stack
    %describes a rectangular pore stack.
    properties (AbortSet=true)
        %% independent properties
        Y1(1,1) double{mustBeNonnegative} %half length of one pore side [m]
        Y2(1,1) double{mustBeNonnegative} %half length of other side    [m]
        Gas_Area_Ratio ...  % percentege of area covered by gas
            (1,1) double{mustBeNonnegative...
            ,mustBeLessThan(Gas_Area_Ratio,1.0001)}=0.1 
        %% dependent geometrical properties
        Diameter(1,1) double{mustBeNonnegative} % stack diameter        [m]
        Area=0.1            % stack area [m^2]         [m]
        A_Solid=0           % Area covered by solid    [m]
        Lp(1,1) double{mustBeNonnegative} %half thickness of solid      [m]
    end
    
    methods
        %% constructor
        function obj = TA_Rectangular_Stack(name,system,length,D,y1,y2,phi,solid)
            obj@TA_Stack(name,system,length,solid)
            obj.Diameter=D;
            obj.Y1=y1;
            obj.Y2=y2;
            obj.Gas_Area_Ratio=phi;
            obj.Lp=(y1*y2*(1/phi-1))/(y1+y2);
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
            obj.Gas_Area_Ratio=1-value/obj.Area;
            
        end
        %% set method for Gas_Area
        function  set.Gas_Area_Ratio(obj,value)
            if value>=1
                error('solid area must be smaller than total area')
            end
            obj.Gas_Area_Ratio=value;
            obj.A_Solid=(1-value)*obj.Area;
        end
        
        %% geometrical functions
        %this function  assists in averaging a property based
        %on the thermal penetration depth and geometryx
        %source:delta E.C version 6.4b2.7 user's guide, equations
        %10.59-10.60
        
        function  f=f_function(obj,delta)
            a=obj.Y1;
            b=obj.Y2;
            %initialize the sum to be added to f
            sum=0;
            for m=1:2:31
                for n=1:2:31
                    YYY=1-1i*pi^2*delta^2/(8*a^2*b^2)*(b^2*m^2+a^2*n^2);
                    sum=sum+1/(m^2*n^2*YYY);    
                end
                f=1-sum*64/pi^4;
            end
        end
        function epsilon=eps_function(obj,deltak,deltas,fk)
            a=obj.Y1;
            b=obj.Y2;
            l=obj.Lp;
            epsilon=fk*(1+1i)*a*b/(deltak*(a+b))/tanh((1+1i)*l/deltas);
            
        end
    end
end

