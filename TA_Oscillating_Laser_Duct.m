
classdef TA_Oscillating_Laser_Duct <TA_Component
    %describes a duct in a thermoacoustic system
    properties (SetObservable=true)
        %% independent properties 
        Length                 % Duct length [m]
        N_Sol_Points           % number of points in solution
        Solid  (1,:)string    % string containing solid name
        A_Solid   (1,1) double ...           %duct solid Area       [m]
            {mustBeNonnegative,mustBeLessThan(A_Solid,1)}
        Q_Laser (1,1) double                    % laser amplitude [W]
        Laserangle (1,1) double                 %laser phase
        D_Beam  (1,1) double                   %laser beam Diameter [m]
        Tout    (1,1) double                   %outside temperature [k]
    end
    properties (Abstract)
        %% dependent geometrical properties
        Area (1,1) double ... % duct area            [m^2]
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}    
        Rh   (1,1) double ...% duct hydraulic radius [m]
            {mustBeGreaterThan(Rh,1e-10),mustBeLessThan(Rh,1)}
                       %heat input
    end
    
    methods
        %% constructor
        function obj = TA_Oscillating_Laser_Duct(name,system,length,solid,q,qangle,Dbeam,Tout)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Solid=solid;            
            obj.Q_Laser=q;
            obj.Laserangle=qangle;
            obj.D_Beam=Dbeam;
            obj.Tout=Tout;
        end
   
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            
            %ode domain
            if isempty(obj.N_Sol_Points)
                domain=[locstart,locstart+obj.Length];
            else
                domain=locstart:obj.Length/(obj.N_Sol_Points-1):locstart+obj.Length;
            end
            ks = solid_properties(obj.Solid,obj.Tout);
            qmlaser=obj.Q_Laser/2;
            Tm=obj.Tout+qmlaser/(2*obj.D_Beam*ks);
            P_U_T_in(3)=Tm;

            [L,P_U_T] = ode45(@(L,P_U_T) obj.ODEfun(L,P_U_T,H_before),domain,P_U_T_in);
            obj.X=L;
            obj.Pressure=P_U_T(:,1);
            obj.Velocity=P_U_T(:,2);
            obj.Temperature=real(P_U_T(:,3));
            PUT_After=P_U_T(end,:);
            PUT_After(3)=real(PUT_After(3));
            H_after=0.5*real(P_U_T(end,1)...
                *conj(P_U_T(end,2)));
        end
        %% ODE method
        function output=ODEfun(obj,~,input,~)
            %this function integrates a radiation driven duct 2020
            
            %input parameter- a vetor containing the complex pressure and
            %velocity and the mean temperature at a point
            %input(1)--p1,input(2)--U1,input(3)--T_m
            %output paremeter-the values of dp/dx, dT/dx and dU/dx at this
            %point
             if ~obj.System_H.Dry_Switch
                error(['Radiation ducts can only',...
                     ' be added to dry systems at the moment'])
             end

            % collect variables from system
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            Array=obj.System_H.Mixture_Array;
            %collect variables from duct
            rh=obj.Rh;
            %A_solid=obj.A_Solid;
            A_gas=obj.Area-obj.A_Solid;
            solid =obj.Solid;
            %collect inputs
            p1=input(1);
            U1=input(2);
            Tm=real(input(3));
            %caclulate mixture and solid properties 
            [ nu, alpha, ~, rho, ~, ~, gamma, ~, ~,cp,k]  = obj.System_H.Mixture_Properties_ODE(p_m,Tm,1,Array);
            deltan=(2*nu/omega)^0.5;
            ks = solid_properties(solid,obj.Tout);
            qm=obj.Q_Laser/2;
            q1=qm*exp(obj.Laserangle);
            %disp(angle(U1))
            Dtube=obj.Rh*4;
            T1=((8*alpha*q1)/(2*obj.D_Beam*ks*pi*Dtube^2)+(1i*omega*p1)/(rho*cp))/(1i*omega+8*alpha/pi/Dtube^2);
            T1s=(q1+2*obj.D_Beam*k*T1)/(2*obj.D_Beam*(ks+k));

            %calculate f function based on different cases
            % based on DEC 10.6-10.7
            if rh/deltan<12.5              %small enough for bessel
                param=(1i-1).*2*rh./deltan; %parameter for convenicence
                
                fnu=(2*besselj(1,param)./(besselj(0,param).*param));            
            elseif rh/deltan>15            %boundary layer
                fnu=(1-1i)*deltan/2/rh;
            else                            %intermediate-interpolation
                param=(1i-1).*2*rh./deltan; 
                fnu1=2*besselj(1,param)./(besselj(0,param).*param);
                fnu2=(1-1i)*deltan/2/rh;
                fnu=fnu1+(rh/deltan-12.5)/(15-12.5)*(fnu2-fnu1);
            end
            
            %

            %F(1)--dP1/dx,F(2)--dU1/dx,F(3)=dTm/dx
             F(1)=-(1i*omega/A_gas*rho/(1-fnu))*U1;
             F(2)=-1i*omega*A_gas/p_m*p1/gamma+8*alpha*(T1s-T1)/(Tm*pi*Dtube^2);
             F(3)=0;

            output=F(:);           
        end
        %% calculate derived properties
        %this function calculates the derived properties (mass flux,
        %acoustic power, total power) based on the pressure, temperature
        % and additional parameters
        function Calculate_Derived(obj)
            if isempty(obj.Pressure)
                error('system must be run before data is collected')
            end
            obj.Acoustic_Power=0.5*real(obj.Pressure.*conj(obj.Velocity));
            
            %calculate total power
            
            if obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
                obj.Total_Power=TP*ones(size(obj.Pressure));
            else
                obj.Total_Power=(obj.Acoustic_Power(1))*ones(size(obj.Pressure));
            end           
            obj.Mass_Flux=0*obj.X;
            obj.Total_Power(end)=obj.Acoustic_Power(end);

            
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a duct, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
        %% get H
        function [H,T1]=get_H(obj)
           error('undefined')
        end
    end
end