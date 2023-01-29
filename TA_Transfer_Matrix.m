classdef  TA_Transfer_Matrix<TA_Component
    % describes an impedance relation between the pressure and velocity
    %before and after. could be used to model a check valve
    %% constant properties
    properties
        Length=0       % length                                 [m]
        N_Sol_Points=2 % number of points in solution
        %% unprotected properties
        M(2,2) double  % transfer matrix
        Jacket...      % indicates whether the surrounding area is jacketed
            (1,1) logical=0;
    end
    
    methods
        function obj=TA_Transfer_Matrix(name,system,M)
            obj@TA_Component(name,system)
            obj.M=M;
        end
          %% set method- to keep length and nsolpoint constant
        function set.Length(~,~)
            error('Length of a Component must be zero')
        end
         function set.N_Sol_Points(~,~)
            error('N_Sol_Points for Component must be 2')
        end
        %% main method
        function [PUT_After,H_after]=run_component(obj,P_U_T_in,H_before,locstart)
            % this function runs the component, getting the variables at the
            %begining as an input, and the variables at the end as an output
            %source-delta E.C version 6.4b2.7 user's guide, section 10.3.1
            %collect variables from system
            PU=P_U_T_in(1:2)';
            PUT_After=[((obj.M+[1,0;0,1])*PU)',P_U_T_in(3)];
            obj.X=locstart*ones(obj.N_Sol_Points,1);
            obj.Pressure=[P_U_T_in(1);PUT_After(1)*ones(1,obj.N_Sol_Points-1)];
            obj.Velocity=[P_U_T_in(2);PUT_After(2)*ones(1,obj.N_Sol_Points-1)];
            obj.Temperature=P_U_T_in(3)*ones(obj.N_Sol_Points,1);
            if obj.Jacket
            H_after=0.5*real(p1_out*conj(U1_out));
            else
                H_after=H_before;
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
            elseif obj.Location>1
                TP=obj.System_H.Components_H{obj.Location-1}.Total_Power(end);
                obj.Total_Power=TP*ones(size(obj.Pressure));
            else
                obj.Total_Power=obj.Acoustic_Power(1)*ones(size(obj.Pressure));
            end
            obj.Mass_Flux=0*obj.X;   
        end
        %% get gas area
        function gas_area=get_gas_area(obj)
            gas_area=NaN;
        end
    end
    
end

