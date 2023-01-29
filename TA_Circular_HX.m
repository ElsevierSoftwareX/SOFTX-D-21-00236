classdef TA_Circular_HX<TA_HX
    %Describes a HX with circular pores
    properties (AbortSet=true)
        %% dependent geometrical properties
        Area=0.1                             % HX area              [m^2]
        Diameter=0.1                         %  Diameter     [m]   
        A_Solid=0;                           % Area covered by solid [m^2]
        Gas_Area_Ratio=0.9;
    end
   
    methods
        %% constructor
        function obj = TA_Circular_HX(name,system,length,D,A,A_solid,solid,heat)
            obj@TA_HX(name,system,length,solid,heat)
            obj.Diameter=D;
            obj.Area=A;
            obj.A_Solid=A_solid;
        end
        
        %% set method for area
        function set.Area(obj,value)
            obj.Area=value;
            obj.A_Solid=value*(1-obj.Gas_Area_Ratio);
        end
        %% set method for A_solid
        function  set.A_Solid(obj,value)
            if value>=obj.Area/2
                error('walls too thick')
            end
            obj.A_Solid=value;
            obj.Gas_Area_Ratio=1-value/obj.Area;
        end
   %% set method for Porosity
        function  set.Gas_Area_Ratio(obj,value)
            if value<=0.5
                error('walls too thick')
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
        function epsilon=eps_function(obj,~,deltas,~)
            rh=obj.Diameter/4;
            epsilon=1/tanh((1+1i)*rh/deltas);
        end
        %% solid temperature calculation
        function T_solid=Calculate_Solid_Temperature(obj)
            % this function calculates the solid temperature based on
            % DEC-delta E.C version 6.4b2.7 user's guide, section 10.7.1
            
            %collect properties from object and system
            dryswitch=obj.System_H.Dry_Switch;
            p_m=obj.System_H.P_m;
            T_m=obj.Temperature(end);
            rh=obj.Diameter/4;
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