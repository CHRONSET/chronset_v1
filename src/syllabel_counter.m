function [nsyl] = syllabel_counter(xidx_p,xidx_n,t,syl_t)

%%
if length(xidx_p) >=1
    d = zeros(1,length(xidx_p));
    k = 0;
    for it = 1:length(xidx_p)
        k = k+1;
        d(1,k) =abs(t(xidx_p(it))-t(xidx_n(it)));
    end;
    d(1,:) = d(1,:)>=syl_t;%minimum syllable duration in time
    
    nsyl= find(d == 1);
    
    if length(nsyl) > length(xidx_p)
        error('number of sylables exceeds number of speech segments');
    end;
    if length(nsyl) > length(xidx_n)
        error('number of sylables exceeds number of speech segments');
    end;
else
    nsyl = [];
end;

