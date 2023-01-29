classdef TA_Round_HX_M<TA_HX
    %Describes a round duct with heat leakege and a water jacket.
    %essentially equivalent to an unideal water jacket
    properties (AbortSet=true)
        %% dependent geometrical properties
        Diameter...                            % HX diameter          [m]
            (1,1) double{mustBeNonnegative}=0.1
        Area=0.1                               % HX area              [m^2]
        Rh=0.1                                 % hydraulic radius     [m]
        Wall_Thickness...                      % thickness of wall    [m]
            (1,1) double{mustBeNonnegative}=0      
        A_Solid=0;                           % Area covered by solid [m^2]
    end
    
    
    methods
        %% constructor
        function obj = TA_Round_HX_M(name,system,length,Diameter,solid,Wall_thickness,heat)
            obj@TA_HX(name,system,length,solid,heat)
            obj.Diameter=Diameter;
            obj.Area=pi*Diameter^2/4;
            obj.A_Solid=pi*Diameter*Wall_thickness;
        end
        %% set method for diameter
        function set.Diameter(obj,value)
            obj.Diameter=value;
            obj.Area=pi*value^2/4;
            obj.Rh=value/4;
            obj.A_Solid=pi*value*obj.Wall_Thickness;
        end
        %% set method for area
        function set.Area(obj,value)
            obj.Area=value;
            %obj.Diameter=eval((value/sym(pi))^0.5*2);
            obj.Diameter=(value/pi)^0.5*2;
            obj.A_Solid=pi*obj.Diameter*obj.Wall_Thickness;
            obj.Rh=obj.Diameter/4;
        end
        %% set method for A_solid
        function  set.A_Solid(obj,value)
            if value>=obj.Area/2
                error('walls too thiick')
            end
            obj.A_Solid=value;
            obj.Wall_Thickness=value/pi/obj.Diameter;
            
        end
        %% set method for wall thickness
        function  set.Wall_Thickness(obj,value)
            if value>=1
                error('solid area must be smaller than total area')
            end
            obj.Wall_Thickness=value;
            obj.A_Solid=value*pi*obj.Diameter;
        end
        %% set method for rh
        function  set.Rh(obj,value)
            obj.Rh=value;
            obj.Diameter=value*4;
            obj.Area=pi*value^2*4;
            obj.A_Solid=pi*obj.Diameter*obj.Wall_Thickness;
        end
        %% geometrical functions
        %this function  assists in averaging a property based
        %on the thermal penetration depth and geometry
        %source:delta E.C version 6.4b2.7 user's guide, equations
        %10.59-10.60
        function  f=f_function(obj,delta)
            rh=obj.Rh;
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
        function epsilon=eps_function(obj,~,deltas,~)
            rh=obj.Rh;
            epsilon=1/tanh((1+1i)*rh/deltas);
        end
        %% solid temperature calculation
        function T_solid=Calculate_Solid_Temperature(obj)
            % this function calculates the solid temperature based on
            % DEC-delta E.C version 6.4b2.7 user's guide, section 10.7.1
            k=50000;
            N=5;
            h=2.5/1000;
            d=8/1000;
            L=obj.Length+0.01;
            q=obj.Heat_Input;
            T=obj.Temperature(1);
            T_solid=T+q*L/h/d/k/N;
        end
    end
    
end