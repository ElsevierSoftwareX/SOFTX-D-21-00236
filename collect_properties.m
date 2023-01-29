function  M_Array=collect_properties(inert,reactive)
%this function collects data from the table for the inert and reactive
%components, into a cell array. the array can be accessed more conveniently
%and with shorter runtimes then a table.
%parameters are collected at the same order they are retrieved in the
%function "mixture_properties"

%input parameters- names of inert and reactive component.
%list of names can be found in the file "property_tables"

if iscell(reactive)
    %% mixture of two reactives  
    M_Array={collect_properties(inert,reactive{1}),...
        collect_properties(inert,reactive{3});reactive{2},reactive{4}};
elseif iscell(inert)
    %% mixture of two inerts
    M_Array={collect_properties(inert{1},reactive),...
        collect_properties(inert{3},reactive);inert{2},inert{4}};  
else
    %% pure gas
    load('property_tables.mat','Property_Tables')
     M_Array=cell(25,length(inert),length(reactive));
    for i=1:length(reactive)
        for j=1:length(inert)
        M_Array{1,j,i}=(Property_Tables.LH{reactive(i),1});
        M_Array{2,j,i}=Property_Tables.Parameters{inert,{'Mw'}};
        M_Array{3,j,i}=Property_Tables.Parameters{reactive(i),{'Mw'}};
        M_Array{4,j,i}=Property_Tables.Parameters{inert,{'miu'}};
        M_Array{5,j,i}=Property_Tables.Parameters{reactive(i),{'miu'}};
        M_Array{6,j,i}=Property_Tables.Parameters{inert,{'Tb'}};
        M_Array{7,j,i}=Property_Tables.Parameters{reactive(i),{'Tb'}};
        M_Array{8,j,i}=Property_Tables.Parameters{inert,{'Tc'}};
        M_Array{9,j,i}=Property_Tables.Parameters{reactive(i),{'Tc'}};
        M_Array{10,j,i}=Property_Tables.Parameters{inert,{'omega'}};
        M_Array{11,j,i}=Property_Tables.Parameters{reactive(i),{'omega'}};
        M_Array{12,j,i}=Property_Tables.Parameters{inert,{'kappa'}};
        M_Array{13,j,i}=Property_Tables.Parameters{reactive(i),{'kappa'}};
        M_Array{14,j,i}=Property_Tables.Parameters{inert,{'a0','a1','a2','a3','a4'}};
        M_Array{15,j,i}=Property_Tables.Parameters{reactive(i),{'a0','a1','a2','a3','a4'}};
        M_Array{16,j,i}=Property_Tables.Parameters{inert,{'Vc'}};
        M_Array{17,j,i}=Property_Tables.Parameters{reactive(i),{'Vc'}};
        M_Array{18,j,i}=Property_Tables.LH{reactive(i),:};
        M_Array{19,j,i}=Property_Tables.TB{reactive(i),:};
        M_Array{20,j,i}=Property_Tables.Parameters{inert,{'ek'}};
        M_Array{21,j,i}=Property_Tables.Parameters{reactive(i),{'ek'}};
        M_Array{22,j,i}=Property_Tables.Parameters{inert,{'sig'}};
        M_Array{23,j,i}=Property_Tables.Parameters{reactive(i),{'sig'}};
        M_Array{24,j,i}=Property_Tables.Parameters{reactive(i),{'beta'}};
        M_Array{25,j,i}=Property_Tables.Parameters{inert,{'beta'}};
        end
    end
end

end

