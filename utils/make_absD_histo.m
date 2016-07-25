function [n,sd] = make_absD_histo(aon,mon,bs)

D = abs(aon-mon);
D = sort(D);

n = zeros(length(bs),1);
sd = zeros(length(bs),1);


for it = 1:length(bs)
    
    if it > 1 && it < length(bs)              
        
        d = [bs(it) bs(it+1)];
        d = abs(diff(d,[],2))/2;
        
        b = [bs(it)-d bs(it+1)-d];
        
        idx = find(D >= b(1) & D <b(2));
        sd(it) = std(D(idx));
        n(it) =length(idx);
    elseif it == 1
        
        d = [bs(it) bs(it+1)];
        d = abs(diff(d,[],2))/2;
        
        b = [bs(it+1)-d];
        
        idx = find(D <b);
        sd(it) = std(D(idx));
        n(it) =length(idx);
    else
        d =(bs(it)-bs(it-1))/2;
        b = bs(it)-d;
        idx = find(D >= b);
        sd(it) = std(D(idx));
        n(it) =length(idx);
        

    end;
    
    
end;


return;