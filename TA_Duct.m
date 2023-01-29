
classdef TA_Duct <TA_Component
    %describes a duct in a thermoacoustic system
    properties (SetObservable=true)
        %% independent properties 
        Length                 % Duct length [m]
        N_Sol_Points           % number of points in solution
        Jacket(1,1) logical=0  % whether or not the duct is water jacketed
    end
    properties (Abstract)
        %% dependent geometrical properties
        Area (1,1) double ... % duct area            [m^2]
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}    
        Rh   (1,1) double ...% duct hydraulic radius [m]
            {mustBeGreaterThan(Rh,1e-10),mustBeLessThan(Rh,1)}  
    end
    
    methods
        %% constructor
        function obj = TA_Duct(name,system,length)
            obj@TA_Component(name,system)
            obj.Length=length;
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
            [L,P_U_T] = ode113A(@(L,P_U_T) obj.ODEfun(P_U_T),domain,P_U_T_in);
            obj.X=L;
            obj.Pressure=P_U_T(:,1);
            obj.Velocity=P_U_T(:,2);
            obj.Temperature=abs(P_U_T(:,3));
            PUT_After=P_U_T(end,:);
            if obj.Jacket
                H_after=0.5*real(P_U_T(end,1)...
                    *conj(P_U_T(end,2)));
            else
                H_after=H_before;
            end
        end
        %% ODE method
        function output=ODEfun(obj,input)
            %this function integrates a duct, including the effects of
            % turbulence, viscocity, thermal and diffusive relaxation
            %the calculation assumes the boundary layer approximation is valid,
            % meaning rh/deltan>30.
            
            %sources:
            % DEC- delta E.C version 6.4b2.7 user's guide, section 10.1.1
            % GST- Gregory swift-Thermoacoustics, A Unifying Perspective,
            %      second edition, chapter 4.4.2. Springer 2017
            % YTP- Theoretical performance characteristics of a travelling-
            %      wave phase-change thermoacoustic engine for low-grade
            %      heat recovery, Rui Yang,Avishai Meir, and Guy Ramon
            %      Applied energy 2020
            
            
            %input parameter- a vetor containing the complex pressure and
            %velocity and the mean temperature at a point
            %input(1)--p1,input(2)--U1,input(3)--T_m
            %output paremeter-the values of dp/dx, dT/dx and dU/dx at this
            %point
            
            % collect variables from system
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            dryswitch=obj.System_H.Dry_Switch;
            Array=obj.System_H.Mixture_Array;
            %collect variables from duct
            A=obj.Area;
            rh=obj.Rh;
            
            %caclulate mixture properties 
            [ nu, alpha, ~, rho, Pr, Sc, gamma, ~, ~,~,~,~,~,~,CB] = obj.System_H.Mixture_Properties_ODE(p_m,input(3),dryswitch,Array);
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
                fnu=fnu1+(rh-12.5)/(15-12.5)*(fnu2-fnu1);
            end
            
            
            %check for turbulence, according to equation DEC 10.15
            NR1=abs(input(2))*rh*4/A/nu;                   %Reynolds number
            if NR1<(500*rh*4/deltan)                       %laminar flow
                m=1;
                m_prime=1;
            else                                            %turbulenct
                [fM,dfMdNR1]=Turbulence_F(5e-4,NR1);       % DEC 10.9
                m=deltan^2*NR1/6/A*(fM-(1-9*pi/32)*NR1*...
                    dfMdNR1)/ (imag(-fnu)/(abs(1-fnu))^2);  % DEC 10.8
                m_prime=((1-deltan/(2*rh))/(1-deltan/(2*rh)/m))^2;  % DEC 10.10
            end
            rnu=omega*rho/A*imag(-fnu)/(abs(1-fnu))^2*m;    % GST 4.74
            ldx=rho/A*( 1-real(fnu))/(abs(1-fnu))^2*m_prime;% GST 4.73
            
            %correct thermal and diffusive relaxation due to turbulence
            deltak=(2*alpha/omega)^0.5*m;
            fk=(1-1i)*deltak/2/rh;
            
            %%% wet term
            if dryswitch~=1
                deltaD=(Pr/Sc)^0.5*deltak;
                fD=(1-1i)*deltaD/2/rh;
                etaD=1;%not precise
                wet1=CB/(1-CB)*fD/etaD*(1-dryswitch); %wet term from YTP 2
            else
                wet1=0;
            end
            % YTP equations 1-3, assuming 0 temperature gradients
            %F(1)--dP1/dx,F(2)--dU1/dx,F(3)=dTm/dx
%             F(1)=-(1i*omega*ldx+rnu)*input(2);

%             F(2)=-1i*omega*A/gamma/p_m*(1+(gamma-1)*fk+gamma*wet1)*input(1);
%             F(3)=0;
            F=[-(1i*omega*ldx+rnu)*input(2);-1i*omega*A/gamma/p_m*(1+(gamma-1)*fk+gamma*wet1)*input(1);0];
            output=F(:);
            %% internal function
            function [fm,dfm] = Turbulence_F(epsilon,Nr)
                %this function returns the value of fm from equation 10.8
                %in the deltaEC user guide. in addition it returns it's derivative with
                %respect to the reynolds number Nr. epsilon is the surface roughness of the
                %pipe, although a value of 5e-4 is always recomended, even for smoother and
                %rougher surfaces
                % the value of f is calculated via equation 10.9 and it's derivative
                %is based on analytic derivation of the equation
                sqrtfm=1;
                for i=1:5
                    sqrtfm=1/(1.7385-2*log10(2*epsilon+18.574/Nr/sqrtfm));
                end
                fm=sqrtfm^2;
                dfm=-32.2663*fm/((1.7385-2*log10(2*epsilon+18.574/(Nr*sqrtfm)))^2*(2*epsilon+18.574/(Nr*sqrtfm))*Nr^2*fm+16.13317*Nr);
            end
            
            
            
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
            if obj.Jacket
                obj.Total_Power=obj.Acoustic_Power;
                if obj.Location>1&&abs(obj.Total_Power(1)-obj.System_H.Components_H{obj.Location-1}.Total_Power(end))>1e-3
                    %issue a warning messege if acoustic power and total
                    %power are different
                    explanation = ['the jacket enforces acoustic power'...
                        ' and total power to be equal. This causes a jump in  '...
                        'acoustic power from the unjacketed area. To prevent this'...
                        ' you can add a TA_Jacket_Target to the system or remove the Jacket.'];
                    hotlinkcode = sprintf('<a href="matlab: disp(''%s'') "> Details </a>', explanation);
                    % Throw warning that includes hyperlink
                    warning('jacket induces jump in acoustic power.  %s', hotlinkcode)
                end
            elseif obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
                obj.Total_Power=TP*ones(size(obj.Pressure));
            else
                obj.Total_Power=obj.Acoustic_Power(1)*ones(size(obj.Pressure));
            end           
            obj.Mass_Flux=0*obj.X;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a duct, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
    end
end

