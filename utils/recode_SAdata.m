function [ID,fn] = recode_SAdata(SAdat)

%%
ID = cell(size(SAdat,1),1);
fn = cell(size(SAdat,1),1);

k = 0;
for it = 1:size(SAdat,1)
        k = k+1;
        ID{k} = SAdat{it,1}(1:regexp(SAdat{it,1},'_')-1);    
        fn{k} = SAdat{it,1}(regexp(SAdat{it,1},'_')+1:regexp(SAdat{it,1},'.W')-1); 
end;

fn = str2double([fn(:)]);
