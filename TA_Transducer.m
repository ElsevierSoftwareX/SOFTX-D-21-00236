classdef  TA_Transducer<TA_Component
    %abstract class, describes an thermoacoustic component
    %% constant properties
    properties
        Length=0            % length                                 [m]
        N_Sol_Points=2      % number of points in solution
    
        %% unprotected properties
        

        Tau (1,1) double   
        Tautag (1,1) double   
        Ze (1,1) double   
        Zm (1,1) double   

        % complex voltage input on loudspeaker   [V]
        V   (1,1) double    
        
        %electric consumption                    [W]
        Electric_Consumption (1,1) double
        
        Jacket  (1,1) logical   %whether or not the loudspeaker is Jacketed
        
    end
    properties (Hidden)
        % an assisting property for using insertsmooth and removesmooth
        %should be removed once links are available
        %when equal to 1 the loudspeaker is "silenced" and has no effect.
        %when equal to 0 the loudspeaker functions regularly.
        Silencer (1,1) double=0;
    end
    methods
        function obj=TA_Transducer(name,system,tau, tautag, ze,zm,v)
            obj@TA_Component(name,system)
            obj.Tau=tau;
            obj.Tautag=tautag;
            obj.Ze=ze;
            obj.Zm=zm;
            obj.V=v;      
        end
          %% set method- to keep length and nsolpoint constant
        function set.N_Sol_Points(~,~)
            error('N_Sol_Points for a Speaker must be 2')
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            %source-delta E.C version 6.4b2.7 user's guide, section 10.3.1
            %collect variables from system
            p_before=P_U_T_in(1);
            U1=P_U_T_in(2);

            %collect variables from object
            tau=obj.Tau;
            tautag=obj.Tautag;
            ze=obj.Ze;
            zm=obj.Zm;
            v=obj.V;  
            
            %ilan this is for you

            p1_out

            
            

            obj.X=locstart*ones(obj.N_Sol_Points,1);
            obj.Pressure=[P_U_T_in(1);p1_out*ones(1,obj.N_Sol_Points-1)];
            obj.Velocity=[P_U_T_in(2);U1_out*ones(1,obj.N_Sol_Points-1)];
            obj.Temperature=P_U_T_in(3)*ones(obj.N_Sol_Points,1);
            PUT_After=[p1_out,U1,T];
            if obj.Jacket
            H_after=0.5*real(p1_out*conj(U1_out));
            else
                H_after=H_before+obj.Electric_Consumption;
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
                if obj.Location>1&&obj.Total_Power(1) ~=obj.System_H.Components_H{obj.Location-1}.Total_Power(end)
                    warning('jacket induces jump in acoustic power. consider removing or adding a jacktarget')        
                end
            else 
                if obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
                else
                TP=obj.Acoustic_Power(1);
                end
                obj.Total_Power=[TP;TP+obj.Electric_Consumption];
            end
            obj.Mass_Flux=0*obj.X;   
        end
        %% addition to Empty variables for Speakers
        function Empty_Varibles(obj)
            obj.Empty_Varibles@TA_Component;
            obj.Electric_Consumption=0;
        end
        %% get gas area
        %this function gets the area of the gas inside the component.
        %In this case of a loudspeaker, gas_area = Area.
        function gas_area=get_gas_area(obj)
            gas_area=ones(length(obj.X),1)*obj.Area;
        end
    end
    
end

