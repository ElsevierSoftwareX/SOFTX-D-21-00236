classdef TA_Circular_Stack <TA_Stack
    %describes a rectangular pore stack.
    properties (AbortSet=true)
        %% independent properties
        Gas_Area_Ratio ...  % percentege of area covered by gas
            (1,1) double{mustBeNonnegative...
            ,mustBeLessThan(Gas_Area_Ratio,1.0001)}=0.1
        Diameter(1,1) double{mustBeNonnegative} %  diameter of the pore [m]

        %% dependent geometrical properties
        Area=0.1            % stack area [m^2]         [m]
        A_Solid=0           % Area covered by solid    [m]
        Lp(1,1) double{mustBeNonnegative}=1e-4 %half thickness of solid      [m]
    end
    
    methods
        %% constructor
        function obj = TA_Circular_Stack(name,system,length,A,D,phi,lp,solid)
            obj@TA_Stack(name,system,length,solid)
            obj.Diameter=D;
            obj.Gas_Area_Ratio=phi;
            obj.Lp=lp;
            obj.Area=A;
            obj.A_Solid=A*(1-obj.Gas_Area_Ratio);
        end
        
        %% set method for area
        function set.Area(obj,value)
            obj.Area=value;
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
        %% set method for Gas_Area_Ratio
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
            rh=obj.Diameter/4;
            if rh/delta<12.5              %small enough for bessel
                param=(1i-1).*2*rh./delta; %parameter for convenicence
                
                f=2*besselj(1,param)./(besselj(0,param).*param);
                
            elseif rh/delta>15            %boundary layer
                f=(1-1i)*delta/2/rh;
            else                            %intermediate-interpolation
                param=(1i-1).*2*rh./delta;
                f1=2*besselj(1,param)./(besselj(0,param).*param);
                f2=(1-1i)*delta/2/rh;
                f=f1+(rh-12.5)/(15-12.5)*(f2-f1);
            end
        end
        function epsilon=eps_function(obj,deltak,deltas,fk)
            r0=obj.Diameter/2;
            lp=obj.Lp;
            epsilon=(fk*(1+1i)*r0/2/deltak)/(tanh((1+1i)*lp/deltas));
        end
    end
end

