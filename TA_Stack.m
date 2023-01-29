classdef TA_Stack <TA_Component
    %abstract component, describes a stack in a TA system
    % subclasses:TA_Plate_Stack, TA_Rectangular_stack
    properties
        %% independent properties
        Length              % Stack length            [m]
        Solid  (1,:)string  % string containing solid name
        N_Sol_Points        % number of points in solution
    end
       properties(Abstract)
        %% Abstract properties
        Area (1,1) double ...                % Stack area            [m^2]
            {mustBeGreaterThan(Area,1e-10),mustBeLessThan(Area,1)}  
        
        A_Solid   (1,1) double ...           %Stack solid Area       [m]
            {mustBeNonnegative,mustBeLessThan(A_Solid,1)}  
    end
    
    methods
        %% constructor
        function obj = TA_Stack(name,system,length,solid)
            obj@TA_Component(name,system)
            obj.Length=length;
            obj.Solid=solid;
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            % begining as an input, and the variables at the end as an output
            
            %ode domain
            if isempty(obj.N_Sol_Points)
                domain=[locstart,locstart+obj.Length];
            else
                domain=locstart:obj.Length/(obj.N_Sol_Points-1):locstart+obj.Length;
            end
            [L,P_U_T] = ode113A(@(L,P_U_T) obj.ODEfun(P_U_T,H_before),domain,P_U_T_in);
            obj.X=L;
            obj.Pressure=P_U_T(:,1);
            obj.Velocity=P_U_T(:,2);
            obj.Temperature=P_U_T(:,3);
            PUT_After=P_U_T(end,:);
            H_after=H_before;
        end
        %% ODE function
        function output=ODEfun(obj,input,H_2)
            %this function integrates a Stack, including the effects of
            % viscocity and the thermoacoustic effect
            
            %source
            %   Theoretical performance characteristics of a travelling-
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
            
            %collect variables from stack
            A=obj.Area;
            A_solid=obj.A_Solid;
            solid =obj.Solid;
            
            %calculate properties -mixture and solid
            [ nu, alpha, D, rho, Pr, Sc, gamma, ~, cp_mol,cp,k,sp,lh,~,CB] = Mixture_Properties(p_m,input(3),dryswitch,Array);
            [k_solid,c_solid,rho_solid] = solid_properties(solid,abs(input(3)));
            beta=1/abs(input(3));%thermal expansion coefficient =1/T
            
            deltanu=(2*nu/omega)^0.5;
            deltak=(2*alpha/omega)^0.5;
            fnu=obj.f_function(deltanu);
            fk=obj.f_function(deltak);
            
            %fnu=tanh((1+1i)*rh/deltanu)/((1+1i)*rh/deltanu);
            %fk=tanh((1+1i)*rh/deltak)/((1+1i)*rh/deltak);
            
            
            delta_s=(2*k_solid/omega/rho_solid/c_solid)^0.5;
            epsilon_s=(k*rho*cp/k_solid/rho_solid/c_solid)^0.5*eps_function(obj,deltak,delta_s,fk);
            Rg=8.314;%J/mol.K
            
            %dry terms
            dry1=1+(gamma-1)*fk/(1+epsilon_s);
            dry2=beta*(fk-fnu)/(1-fnu)/(1-Pr)/(1+epsilon_s);
            dry3=0.5*real(input(1)*conj(input(2))*(1-(fk-conj(fnu))/(1+epsilon_s)/(1+Pr)/(1-conj(fnu))));
            dry4=-cp_mol/(1-Pr)*imag(conj(fnu)+(fk-conj(fnu))*(1+epsilon_s*fnu/fk)/(1+epsilon_s)/(1+Pr));

            if dryswitch~=1
                deltaD=(2*D/omega)^0.5;
                fD=obj.f_function(deltaD);
 
                etaD=1;%not precise
                etanu=1;%not precise
                wet1=CB/(1-CB)*fD/etaD*(1-dryswitch);
                wet2=(fnu-etanu/etaD*fD)*CB*lh/(1-fnu)/(1-Sc)/(1-CB)/Rg/abs(input(3))^2*(1-dryswitch);
                wet3=CB/(1-CB)*lh/Rg/abs(input(3))/(1+Sc)*real(conj(input(2))*input(1)/(A-A_solid)*(1-fD-(1-conj(fnu)))/(1-conj(fnu)))*(1-dryswitch);
                wet4=CB/(1-CB)*lh^2/Rg/abs(input(3))^2/(1-Sc^2)*imag((1-conj(fnu))*(1+Sc)-etanu/etaD*(1+fD+1-conj(fnu)))*(1-dryswitch);
                wet5=D*lh^2*CB/(1-CB)/Rg/abs(input(3))^2*(1-dryswitch);             
            else
                wet1=0;
                wet2=0;
                wet3=0;
                wet4=0;
                wet5=0;
            end
            %F(1)--dP1/dx,F(2)--dU1/dx,F(3)=dTm/dx
            %equations 1-3
            F(1)=-1i*omega*rho*input(2)/(1-fnu)/(A-A_solid);
            F(3)=(-H_2+dry3+0.5*wet3*(A-A_solid)) / (p_m/Rg/abs(input(3))*((abs(input(2)))^2/2/(A-A_solid)/omega/(abs(1-fnu))^2*(dry4+wet4)+wet5*(A-A_solid))+((A-A_solid)*k+A_solid*k_solid));
            F(2)=-1i*omega*(A-A_solid)/rho/sp^2*(dry1+gamma*wet1)*input(1)+(dry2-wet2)*F(3)*input(2);%
            
            output=F(:);
        end
        %% calculate derived properties
        function Calculate_Derived(obj)
            if isempty(obj.Pressure)
                error('system must be run before data is collected')
            end
            obj.Acoustic_Power=0.5*real(obj.Pressure.*conj(obj.Velocity));
            if obj.Location>1
                obj.Total_Power=obj.System_H.Components_H{obj.Location-1}.Total_Power(end)*ones(length(obj.Acoustic_Power),1);
            else
                obj.Total_Power=obj.Acoustic_Power(1)*ones(size(obj.Pressure));
            end
            %caculate mass flux
            obj.Mass_Flux=0*obj.Pressure;
            for i=1:length(obj.Pressure)
                obj.Mass_Flux(i,1)=Get_Mass_Flux(obj.Pressure(i),obj.Velocity(i),obj.Temperature(i));
            end
            
            %internal function
            function massflux=Get_Mass_Flux(p1,U1,T)
                %this function calculates the mass flux based on:
                %Theoretical performance characteristics of a travelling-
                %      wave phase-change thermoacoustic engine for low-grade
                %      heat recovery, Rui Yang,Avishai Meir, and Guy Ramon
                %      Applied energy 2020, equation 6
                dryswitch=obj.System_H.Dry_Switch;
                if dryswitch
                    massflux=0;
                    return
                end
                % collect variables from system
                omega=obj.System_H.Frequency*2*pi;
                p_m=obj.System_H.P_m;
                Array=obj.System_H.Mixture_Array;
                solid =obj.Solid;
                
                %collect variables from stack
                A=obj.Area;
                A_solid=obj.A_Solid;
                H_2=obj.Total_Power(1);
                %calculate mixture and solidproperties
                [ nu, kappa, D, rho, ~, Sc, ~, ~, cp_mol,cp,k,~,lh,~,CB,lh_m] = Mixture_Properties(p_m,T,dryswitch,Array);
                [k_solid,c_solid,rho_solid] = solid_properties(solid,T);
                
                deltanu=(2*nu/omega)^0.5;
                deltak=(2*kappa/omega)^0.5;
                deltaD=(2*D/omega)^0.5;
                fnu=obj.f_function(deltanu);
                fk=obj.f_function(deltak);
                fD=obj.f_function(deltaD);
                
                
                Pr=(deltanu/deltak)^2;
                delta_s=(2*k_solid/omega/rho_solid/c_solid)^0.5;
                epsilon_s=(k*rho*cp/k_solid/rho_solid/c_solid)^0.5*obj.eps_function(deltak,delta_s,fk);
                
                Rg=8.314;%J/mol.K
                
                
                etaD=1;%not precise
                etanu=1;%not precise
                
                wet3=CB/(1-CB)*lh/Rg/abs(T)/(1+Sc)*real(conj(U1)*p1/(A-A_solid)*(1-fD-(1-conj(fnu)))/(1-conj(fnu)));
                dry3=0.5*real(p1*conj(U1)*(1-(fk-conj(fnu))/(1+epsilon_s)/(1+Pr)/(1-conj(fnu))));
                wet4=CB/(1-CB)*lh^2/Rg/abs(T)^2/(1-Sc^2)*imag((1-conj(fnu))*(1+Sc)-etanu/etaD*(-1+fD+1-conj(fnu)));
                dry4=-cp_mol/(1-Pr)*imag(conj(fnu)+(fk-conj(fnu))*(1+epsilon_s*fnu/fk)/(1+epsilon_s)/(1+Pr));
                wet5=D*lh^2*CB/(1-CB)/Rg/abs(T)^2;
                
                dtdx=(-H_2+dry3+0.5*wet3*(A-A_solid)) / (p_m/Rg/abs(T)*((abs(U1))^2/2/(A-A_solid)/omega/(abs(1-fnu))^2*(dry4+wet4)+wet5*(A-A_solid))+((A-A_solid)*k+A_solid*k_solid));
                
                heatflux=0.5*wet3*(A-A_solid)-dtdx*(p_m/Rg/abs(T)*((abs(U1))^2/2/(A-A_solid)/omega/(abs(1-fnu))^2*(wet4)+wet5*(A-A_solid))); %j/s
                massflux=heatflux/lh_m; %kg/s
%                 heatflux-0.5*wet3*(A-A_solid);
%                 dtdx*(CB/(1-CB)*lh^2/Rg/abs(T)^2/(1-Sc^2)*p_m/Rg/abs(T)*((abs(U1))^2/2/(A-A_solid)/omega/(abs(1-fnu))^2 ...
%                     *(imag((1-conj(fnu)) ...
%                     *(1+Sc)+(1-fD-1+conj(fnu)))))); %j/s
            end
        end
        %% analyze performance
        function result=Analyze(obj)
            %this function returns data regarding the stack's performance
            if isempty(obj.Acoustic_Power)
                error('run system while collecting data before analyzing')
            end
            [term1,term2,term3]=deal(0*obj.Pressure);
            omega=obj.System_H.Frequency*2*pi;
            p_m=obj.System_H.P_m;
            dryswitch=obj.System_H.Dry_Switch;
            Array=obj.System_H.Mixture_Array;
            Rg=8.314;%J/mol.K
            dx=diff(obj.X);
            if max(dx)-min(dx)<1e-10
                dx=dx(1);
            else 
                error('cannot derive temperature in variable X stack')
            end
            f=obj.Temperature;
            DTDX=0*f;
            DTDX(2:end-1)=(f(3:end)-f(1:end-2))/(2*dx);
            DTDX(1)=(-1.5*f(1)+2*f(2)-0.5*f(3))/dx;
            DTDX(end)=(1.5*f(end)-2*f(end-1)+0.5*f(end-2))/dx;
            for i=1:length(obj.Pressure)
               dtdx=DTDX(i);
               T_m=obj.Temperature(i);
               [ nu, ~, D, ~, ~, Sc, ~, ~, ~,~,~,~,lh,...
                   ~,CB,lh_m] = Mixture_Properties(p_m,T_m,dryswitch,Array);
                deltanu=(2*nu/omega)^0.5;
                deltaD=(2*D/omega)^0.5;
                p1=obj.Pressure(i);
                U1=obj.Velocity(i);
                Fnu=1-obj.f_function(deltanu);
                FD=1-obj.f_function(deltaD);
                term1(i)=CB/(1-CB)/2/Rg/T_m*real(conj(U1)*p1*(FD-conj(Fnu))...
                    /(conj(Fnu)*(1+Sc)))*lh;
                term2(i)=1/(1-CB)*abs(U1)^2/2/obj.Area/obj.Gas_Area_Ratio/omega/abs(Fnu)^2 ...
                *p_m*CB*lh/Rg^2/T_m^3*1/(1-Sc^2)*imag(conj(Fnu)*(1+Sc)+FD-conj(Fnu))*dtdx*lh;
                term3(i)=p_m*obj.Area*obj.Gas_Area_Ratio*D/(1-CB)*CB*lh_m/Rg^2/T_m^3*dtdx*lh;
            end
            result.Mass_Flux_Terms=[term1,term2,term3];
            result.Impedance=abs(obj.Pressure)./abs(obj.Velocity)*obj.Area*obj.Gas_Area_Ratio;
            result.Temperature_Gradient=DTDX;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a stack, gas_area = Area - Solid_Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*(obj.Area-obj.A_Solid);
        end
    end
    methods(Abstract)
        %each stack must have functions to represent the geometry of the
        % stack, in term of the parameter f and epsilon.
        % see-delta E.C version 6.4b2.7 user's guide, section 10.5.1
        
        %f function, representing spactial average inside pore
        f=f_function(obj,delta)
        %epsilon function, representing solid activity divided by the fluid
        %(always identical) first term.
        e=eps_function(obj,deltak,deltas,fk)
    end
    
end

