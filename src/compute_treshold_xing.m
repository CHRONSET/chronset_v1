function [x_tresh,x_idx] = compute_treshold_xing(features,tresh)

%%
x_tresh = cell(1,length(features));
x_idx = cell(1,length(features));
%%
for it = 1:length(features)
    
    if ismember(it,[1 3 4 6]) 
        
        x_idx = find(features{it}>tresh{it});        
        x_tresh{it} = zeros(1,length(features{it}));
        x_tresh{it}(x_idx) =1;
        
    elseif ismember(it,[2 5]) 
        
        x_idx = find(features{it}<tresh{it});        
        x_tresh{it} = zeros(1,length(features{it}));
        x_tresh{it}(x_idx) =1;
        
    end;    
    
end;