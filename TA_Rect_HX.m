classdef TA_Rect_HX <TA_HX
    %Describes a Rectanbular
    properties(AbortSet=true)
        %% independent properties
        Y1(1,1) double{mustBeNonnegative} %half length of one pore side [m]
        Y2(1,1) double{mustBeNonnegative} %half length of other side    [m]
        Gas_Area_Ratio ...  % percentege of area covered by gas
            (1,1) double{mustBeGreaterThan(Gas_Area_Ratio,0)...
            ,mustBeLessThan(Gas_Area_Ratio,1)}  =0.1      
         %% dependent geometrical properties
        Diameter (1,1) double{mustBeNonnegative}  % HX diameter         [m]
        Area=0.1                % HX area [m^2]         [m]
        A_Solid =0.1            % Area covered by solid    [m]
        Lp(1,1) double{mustBeNonnegative} %half thickness of solid      [m]
    end
    
    methods
        %% constructor
        function obj = TA_Rect_HX(name,system,length,D,y1,y2,phi,solid,Q)
            obj@TA_HX(name,system,length,solid,Q)
            obj.Diameter=D;
            obj.Y1=y1;
            obj.Y2=y2;
            obj.Gas_Area_Ratio=phi;
            obj.Lp=(y1*y2*(1/phi-1))/(y1+y2);
            obj.Area=pi*D^2/4;
            obj.A_Solid=pi*D^2/4*obj.Gas_Area_Ratio;
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
        %% solid temperature calculation
        function T_solid=Calculate_Solid_Temperature(obj)
            % this function calculates the solid temperature based on
            % DEC-delta E.C version 6.4b2.7 user's guide, section 10.7.1
            
            %collect properties from obje ct and system
            dryswitch=obj.System_H.Dry_Switch;
            p_m=obj.System_H.P_m;
            T_m=obj.Temperature(end);
            rh=(obj.Y1*obj.Y2)/2/(obj.Y1+obj.Y2);
            Array=obj.System_H.Mixture_Array;
            omega=obj.System_H.Frequency*2*pi;
            
            %calculate properties
            [ ~, alpha, ~, ~, ~, ~, ~, ~, ~,~,k] = Mixture_Properties(p_m,T_m,dryswitch,Array);
            
            zeta= abs(obj.Velocity(end)/(obj.Area-obj.A_Solid))/omega;              %acoustic displacement
            xeff=min(2*zeta,obj.Length);
            
            deltak=(2*alpha/omega)^0.5;
            yeff=min(deltak,rh);                                                %effective distance from the plate
            h=k/yeff;                                                           %approximate convection coefficient
            T_solid=T_m+obj.Heat_Input/(h*(obj.Area-obj.A_Solid)*xeff/rh);      %DEC eq 10.113
            
        end
    end
end

