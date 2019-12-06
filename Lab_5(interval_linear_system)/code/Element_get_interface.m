function [result] = Element_get_interface()
    result = struct('create', @Element_Create, 'print', @Element_Print, 'draw', @Element_Draw, 'get_cut', @Element_GetCut, 'reindexing', @Element_Reindexing, 'get_hord', @Element_GetHord);
    
end