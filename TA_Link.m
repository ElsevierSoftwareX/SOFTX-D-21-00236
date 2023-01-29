classdef TA_Link<matlab.mixin.Copyable
%this class defines a link between two or more properties in a TA_System    
    properties (SetAccess=?TA_System)
        Properties (1,:)cell
    end
     properties (SetAccess=?TA_System)
        System (1,:) char  
        Name   (1,:) char
    end
    properties (Hidden,SetAccess=?TA_System)
       System_H TA_System
    end
    properties (SetAccess=?TA_System)
       Listeners (1,:) event.proplistener
    end

    properties 
       Linked_Properties (1,:) cell % what is this definition? can it define infinite table?
       %Listeners
    end
   
    methods
    %% constructor 

    function obj = TA_Link(system,Name)
            obj.System = system.Name;
            obj.Name=Name;
            obj.System_H = system;

            obj.Linked_Properties = {};

        end

        function Add_Link(obj,Linked_Property,Linked_Components_List)
            obj.Linked_Properties{1} = Linked_Property;
            for i = 1:length(Linked_Components_List)
                obj.Listeners(i) = addlistener(obj.System_H.Components_H{i},Linked_Property,'PostSet',@obj.Equal_Linked_Properties); 
                obj.Linked_Properties{i+1} = Linked_Components_List(i);
            end
        end
        
        function Equal_Linked_Properties(obj,Source,Data)
            
            for i = 1:length(obj.Linked_Properties)-1
                
                %obj.System_H.Components_H{i}.(obj.Linked_Properties{1})

                obj.Listeners(i).Enabled = false;
                obj.System_H.Components_H{i}.(obj.Linked_Properties{1}) = Data.AffectedObject.(Data.Source.Name);
                obj.Listeners(i).Enabled = true;
            end
           
        end

        function Delete_Link(obj,Linked_Property,Delete_Components_Link)
            for i = 1:length(Delete_Components_Link)
                delete(obj.Listeners(Delete_Components_Link(i)))
                obj.Listeners(Delete_Components_Link(i)) = [];
            end
        end
    end
end

