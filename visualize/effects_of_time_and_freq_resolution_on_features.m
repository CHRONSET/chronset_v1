%%
figure;
g = 1 - (size(savedata,2)/(size(savedata,2)*2));

for it =1:size(savedata,1)
    
    k=0;
    for kt = 1:length(savedata{it,1}.features)
        k = k+1;
        c = [0 0 g];
        for jt=1:size(savedata,2)
            c = c+.1;
            subplot(3,2,k);
            hold on;
            plot(savedata{it,jt}.features{kt},'Color',c);
            
        end;
    end;
end;