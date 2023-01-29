
classdef TA_Radiation_Duct <TA_Component
    %describes a duct in a thermoacoustic system
    properties (SetObservable=true)
        %% independent properties 
        Length                 % Duct length [m]
        N_Sol_Points           % number of points in solution
        Solid  (1,:)string    % string containing solid name
        A_Solid   (1,1) double ...           %duct solid Area       [m]
            {mustBeNonnegative,mustBeLessThan(A_Solid,1)}
        Absorptivity (1,1) double ...           %duct solid Area       [m]
            {mustBeNonnegative}
        Q_in (1,1) double 
        Jacket (1,1) logical=0
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
        function obj = TA_Radiation_Duct(name,system,length,kappa,solid,q)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Absorptivity=kappa;
            obj.Solid=solid;            
            obj.Q_in=q;
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
            [L,P_U_T] = ode45(@(L,P_U_T) obj.ODEfun(L,P_U_T,H_before+obj.Q_in),domain,P_U_T_in);
            obj.X=L;
            obj.Pressure=P_U_T(:,1);
            obj.Velocity=P_U_T(:,2);
            obj.Temperature=real(P_U_T(:,3));
            PUT_After=P_U_T(end,:);
            PUT_After(3)=real(PUT_After(3));
            if obj.Jacket
            H_after=0.5*real(P_U_T(end,1)...
                    *conj(P_U_T(end,2)));
            else
                H_after=H_before+obj.Q_in;
            end
        end
        %% ODE method
        function output=ODEfun(obj,~,input,H)
            %this function integrates a radiation driven duct 2020
            
            %disp(X-0.1)
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
            kappa=obj.Absorptivity;
            sigma=5.670374419e-8; %[W/m^2k^4]
            Array=obj.System_H.Mixture_Array;
            %collect variables from duct
            rh=obj.Rh;
            A_solid=obj.A_Solid;
            A_gas=obj.Area-obj.A_Solid;
            solid =obj.Solid;
            %collect inputs
            p1=input(1);
            U1=input(2);
            Tm=real(input(3));
            %caclulate mixture and solid properties 
            [ nu, ~, ~, rho, ~, ~, ~, ~, ~,cp] = obj.System_H.Mixture_Properties_ODE(p_m,Tm,1,Array);
            [k_solid] = solid_properties(solid,abs(Tm));
            deltan=(2*nu/omega)^0.5;
            
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
     
            %F(1)--dP1/dx,F(2)--dU1/dx,F(3)=dTm/dx
             F(1)=-(1i*omega/A_gas*rho/(1-fnu))*U1;
             F(3)=(H-rho*real(conj(U1)*p1*rho+16i*kappa*sigma*Tm^3*conj(U1)*p1/cp/omega)...
                 /(2*(rho^2+(16*kappa*sigma*Tm^3/cp/omega)^2)))...
                 /((-8*rho^2*kappa*sigma*Tm^3*abs(U1)^2)/...
                 (A_gas*omega^2*(rho^2+(16*kappa*sigma*Tm^3/cp/omega)^2))...
                 -A_solid*k_solid-A_gas*16*sigma*Tm^3/3/kappa);
             F(3)=real(F(3));
             F(2)=-1i*omega*A_gas/p_m*p1+...
                 (1i*omega*A_gas*p1-1i*16*kappa*sigma*Tm^3/omega*F(3)*U1)...
                 /(Tm*rho*cp-1i*16*kappa*sigma*Tm^4/omega);%
            if (F(3))>0
                a=0;
            end
            %disp(imag(Tm))
            %disp(imag(F(3)))
            output=F(:);
            %disp(Tm)
            %disp(U1)
            %disp(F(2))
            if(imag(F(3)))>0||imag(Tm)>0
                nathan=9
            end
            %ratio=abs((-8*rho^2*kappa*sigma*Tm^3*abs(U1)^2)/...
               %  (A_gas*omega^2*(rho^2+(16*kappa*sigma*Tm^3/cp/omega)^2)))...
               %  /abs(-A_solid*k_solid)    ;
            
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
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end)+obj.Q_in;
                obj.Total_Power=TP*ones(size(obj.Pressure));
            else
                obj.Total_Power=(obj.Acoustic_Power(1)+obj.Q_in)*ones(size(obj.Pressure));
            end           
            obj.Mass_Flux=0*obj.X;
            if obj.Jacket
            obj.Total_Power(end)=obj.Acoustic_Power(end);
            end
            
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a duct, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
        %% get H
        function [H,T1]=get_H(obj)
            kappa=obj.Absorptivity;
            sigma=5.670374419e-8; %[W/m^2k^4]
            omega=obj.System_H.Frequency*2*pi;
            A_gas=(obj.Area-obj.A_Solid);
            dTdx=obj.System_H.diff12(obj.Temperature,obj.X);
            T1=zeros(size(obj.Pressure));
            H=T1;
            for i=1:length(T1)
            [ ~, ~, ~, rho, ~, ~, ~, ~ ,~,cp] = Mixture_Properties(obj.System_H.P_m,obj.Temperature(i),obj.System_H.Dry_Switch,obj.System_H.Mixture_Array);
            [ks,~,~]=solid_properties(obj.Solid,obj.Temperature(i));
            p1=obj.Pressure(i);
            u1=obj.Velocity(i)/A_gas;
            T1(i)=(p1/cp+1i*rho*u1/omega*dTdx(i))/(rho-1i*16*kappa*sigma/omega/cp*obj.Temperature(i)^3);
            H(i)=A_gas*rho*cp*0.5*real(T1(i)*conj(u1))-obj.A_Solid*ks*dTdx(i)...
                -A_gas*16/3*obj.Temperature(i)^3/kappa*sigma/rho/cp*dTdx(i);
            end
            %H=0;
        end
    end
end

