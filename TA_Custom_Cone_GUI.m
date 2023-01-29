classdef TA_Custom_Cone_GUI <TA_Cone
    %describes a duct in a thermoacoustic system
    properties (AbortSet=true)
        %% dependent geometrical properties
       
        Area_1      (1,1) double{mustBeNonnegative} % first area     [m^2]
        Area_2      (1,1) double{mustBeNonnegative} % second area    [m^2]
        Rh_1        (1,1) double{mustBeNonnegative}  % hydraulic radiuses [m]
        Rh_2        (1,1) double{mustBeNonnegative}
        ModeR % linear,exponent+, exponent-, square, sqrt
        ModeA % linear,exponent+, exponent-, square, sqrt
    end   
    
    methods
        %% constructor
        function obj = TA_Custom_Cone_GUI(name,system,length,Rh_1,Rh_2,Area_1,Area_2,ModeR,ModeA)
            obj@TA_Cone(name,system,length)
            obj.Area_1=Area_1; 
            obj.Area_2=Area_2;
            obj.Rh_1=Rh_1;
            obj.Rh_2=Rh_2;
            obj.ModeR=ModeR;
            obj.ModeA=ModeA;
        end
        %% calculation methods
        function A=Calc_A(obj,x_rel)
            switch obj.ModeA
                case 'linear' %y=a*x+b
                    A=obj.Area_1+(obj.Area_2-obj.Area_1)*x_rel/obj.Length;
                case 'exponent+' % y=a+b*exp(x)
                    A=obj.Area_1+((obj.Area_2-obj.Area_1)*(1-exp(x_rel)))/(1-exp(obj.Length));
                case 'exponent-'% y=a+b*exp(-x)
                    A=obj.Area_1+((obj.Area_2-obj.Area_1)*(1-exp(-x_rel)))/(1-exp(-obj.Length));
                case 'sqrt' % y=a+bx^0.5
                    A=obj.Area_1+((obj.Area_2-obj.Area_1)/obj.Length^0.5)*x_rel^0.5;
                case 'square' % y=a+bx^2
                    A=obj.Area_1+((obj.Area_2-obj.Area_1)/obj.Length^2)*x_rel^2;
                otherwise 
                    error('wrong input format')
            end
        end
        function rh=Calc_Rh(obj,x_rel)
            switch obj.ModeR
                 case 'linear' %y=a*x+b
                    rh=obj.Rh_1+(obj.Rh_2-obj.Rh_1)*x_rel/obj.Length;
                case 'exponent+' % y=a+b*exp(x)
                    rh=obj.Rh_1+((obj.Rh_2-obj.Rh_1)*(1-exp(x_rel)))/(1-exp(obj.Length));
                case 'exponent-' % y=a+b*exp(-x)
                    rh=obj.Rh_1+((obj.Rh_2-obj.Rh_1)*(1-exp(-x_rel)))/(1-exp(-obj.Length));
                case 'sqrt' % y=a+bx^0.5
                    rh=obj.Rh_1+((obj.Rh_2-obj.Rh_1)/obj.Length^0.5)*x_rel^0.5;
                case 'square' % y=a+bx^2
                    rh=obj.Rh_1+((obj.Rh_2-obj.Rh_1)/obj.Length^2)*x_rel^2;
                otherwise 
                    error('wrong input format')
            end
        end
    end
end

