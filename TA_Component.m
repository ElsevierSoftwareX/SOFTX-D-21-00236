classdef  TA_Component<matlab.mixin.Copyable
    %abstract class, describes an thermoacoustic component
    
    properties (Abstract)
        %% unprotected properties
        % length [m]
        Length (1,1) double{mustBeNonnegative}
        % number of points in solution
        N_Sol_Points double{mustBeInteger}
    end
    properties (SetAccess=?TA_System)
        %% protected properties
        % string with name of "mother" system of the component
        System (1,:) char
        
        % the location of the duct in the mother system
        Location double{mustBeInteger}
    end
    properties(SetAccess= 'protected')
        %% variable profile
        %these are the governing equation's parameters:oscilating pressure,
        %oscillationg velocity, mean temperature Location and some derived
        % properties.
        
        %string containing variable name
        Name           (1,:) char
        
        X              (:,1) double   %the location                   [m]
        Pressure       (:,1) double   %the oscillating pressure field [pa]
        Velocity       (:,1) double   %the oscillating velocity field [m/s]
        Temperature     ...           %temperature                    [k]
            (:,1) double{mustBeNonnegative}
        
        Acoustic_Power (:,1) double   %the acoustic power [W]
        Total_Power    (:,1) double   % the total power   [W]
        Mass_Flux      (:,1) double   % the mass flux     [W]
        
    end
    properties (Hidden,SetAccess=?TA_System)
        %% internal parameters for the class to work
        System_H  TA_System %handle to actual system
    end
    methods
        %% constructor
        function obj=TA_Component(name,system)
            obj.Name=name;
            if ~(isa(system,'char')&&strcmp(system,'null'))
                obj.System=system.Name;
                obj.System_H=system;
                obj.Location=size(system.Components,2)+1;
                system.addComponent(name);
                system.Components_H=[system.Components_H,{obj}];
            end
        end
        %% rename function
        function rename(obj,value)
            oldName=obj.Name;           
            obj.Name=value;
            if ~isempty(obj.System)
                [obj.System_H.Guesses]=obj.System_H. ...
                SearchArray(obj.System_H.Guesses,oldName,value);
                [obj.System_H.Targets]=obj.System_H.SearchArray ...
                (obj.System_H.Targets,oldName,value);
                obj.System_H.Components{obj.Location}=value;
            end
            
            
        end
        %% copy component
        function cp = Copy_Component(obj,varargin)
            %copy a component
            % Copy_Component() creates a copy of the component with the
            %                  name "copy of 'old name'"
            % Copy_Component(newName) creates a copy of the component with
            % the name specified
            
            % Shallow copy object
            cp = obj.copyElement;
            cp.System_H=TA_System.empty;
            cp.System=[];
            cp.Location=[];
            cp.Empty_Varibles;
            
            if isempty(varargin)
                cp.Name=['copy of ',cp.Name];
            else
                cp.Name=varargin{1};
            end
        end
        %% empty all variables (presure, velocity etc)
        function Empty_Varibles(obj)
            obj.X=[];
            obj.Pressure=[];
            obj.Velocity=[];
            obj.Temperature=[];
            obj.Acoustic_Power=[];
            obj.Total_Power=[];
            obj.Mass_Flux=[];
        end
        %% caclulate Fs
        function [FD,Fv,Falpha] = CalculateFs(obj)
            %this function calculates the F functions FD Fv Falpha for a component
            % the F functions are defined based on Offner's notation (JFM 2019), and
            % are equal to 1-f based on the classical thermoacoustic notation
            
            % check for errors
            if isempty(obj.Pressure)
                error('Run System before collecting Fs')
            end
            
            % collect system data
            Sys= obj.System_H;
            M_Array=Sys.Mixture_Array;
            dry_switch=Sys.Dry_Switch;
            P_m=Sys.P_m;
            Omega=Sys.Frequency*2*pi;
            % assign fs
            
            [FD,Fv,Falpha]=deal(ones(length(obj.X),1));        %initialise variable
            
            if isa(obj,'TA_Duct')||isa(obj,'TA_Cone')  % check type of component
                % ducts and cones
                for i=1:length(obj.X)
                    %collect rh
                    if isa(obj,'TA_Duct')
                        rh=obj.Rh;                         %duct hydraulic radius
                    else
                        rh=obj.Calc_Rh(obj.X(i));%cone hydraulic radius
                    end
                    %caclulate mixture properties
                    [ nu, alpha, ~, ~, Pr, Sc] = Mixture_Properties(P_m,obj.Temperature(i),dry_switch,M_Array);
                    deltan=(2*nu/Omega)^0.5;
                    deltak=(2*alpha/Omega)^0.5;
                    deltaD=(Pr/Sc)^0.5*deltak;
                    
                    %calculate f function based on different cases
                    % based on DEC 10.6-10.7
                    if rh/deltan<12.5              %small enough for bessel
                        paramNU=(1i-1).*2*rh./deltan; %parameter for convenicence
                        paramK=(1i-1).*2*rh./deltaK;
                        paramD=(1i-1).*2*rh./deltaD;
                        
                        fnu=(2*besselj(1,paramNU)./(besselj(0,paramNU).*paramNU));
                        fk=(2*besselj(1,paramK)./(besselj(0,paramK).*paramK));
                        fD=(2*besselj(1,paramD)./(besselj(0,paramD).*paramD));
                    elseif rh/deltan>15            %boundary layer
                        fnu=(1-1i)*deltan/2/rh;
                        fk=(1-1i)*deltak/2/rh;
                        fD=(1-1i)*deltaD/2/rh;
                    else                            %intermediate-interpolation
                        paramNU=(1i-1).*2*rh./deltan; %parameter for convenicence
                        paramK=(1i-1).*2*rh./deltaK;
                        paramD=(1i-1).*2*rh./deltaD;
                        
                        fnu1=2*besselj(1,paramNU)./(besselj(0,paramNU).*paramNU);
                        fk1=2*besselj(1,paramK)./(besselj(0,paramK).*paramK);
                        fD1=2*besselj(1,paramD)./(besselj(0,paramD).*paramD);
                        
                        fnu2=(1-1i)*deltan/2/rh;
                        fk2=(1-1i)*deltak/2/rh;
                        fD2=(1-1i)*deltaD/2/rh;
                        
                        fnu=fnu1+(rh-12.5)/(15-12.5)*(fnu2-fnu1);
                        fk=fk1+(rh-12.5)/(15-12.5)*(fk2-fk1);
                        fD=fD1+(rh-12.5)/(15-12.5)*(fD2-fD1);
                    end
                    FD(i)=1-fD;
                    Falpha(i)=1-fk;
                    Fv(i)=1-fnu;
                end
            elseif isa(obj,'TA_Stack')||isa(obj,'TA_HX')
                for i=1:length(obj.X)
                    [ nu, alpha, ~, ~, Pr, Sc] = Mixture_Properties(P_m,obj.Temperature(i),dry_switch,M_Array);
                    deltan=(2*nu/Omega)^0.5;
                    deltak=(2*alpha/Omega)^0.5;
                    deltaD=(Pr/Sc)^0.5*deltak;
                    
                    FD(i)=1-obj.f_function(deltaD);
                    Falpha(i)=1-obj.f_function(deltak);
                    Fv(i)=1-obj.f_function(deltan);
                end
            end
        end
        
        
    end
    
    methods (Abstract)
        run_component(obj,P_U_T_in,H_in,locstart)
        Calculate_Derived(obj)
        get_gas_area(obj)
    end
end

