classdef TA_Equals_Link<TA_Link
%this class defines two or more properties in a TA_System as equal.   
    properties (SetAccess=?TA_System)
        Properties (1,:)cell
    end
    properties (Hidden,SetAccess=?TA_System)
        Listeners (1,:) event.proplistener
    end
    methods
    %% constructor 
        function obj = TA_Equals_Link(system,Name)
         obj@TA_Link(system,Name)
        end

        function (TASource,Taproperty)
            addlistener(TASource,Taproperty,'PostSet',@obj.Equals_Func);
            addlistener(TASource,'Diameter','PostSet',@obj.Equals_Func);
            addlistener(TASource,'DiameterChangedEvent',@obj.ListenerReaction)
        end
        %% equalise
        %change all the 
    end
end

