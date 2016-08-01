function joint_distribution_plot_measurement_errors(d1,d2,D)
%%
x_bin = linspace(min(min(D))-10,max(max(D))+10,round((round(max(max(D))-min(min(D)))*1)/50));
y_bin = linspace(min(min(D))-10,max(max(D))+10,round((round(max(max(D))-min(min(D)))*1)/50));

x_bin = round(x_bin*10)/10;
y_bin = round(y_bin*10)/10;
%%
M = [d1 d2];
[N] = zeros(length(y_bin),length(x_bin));
for it = 1:length(x_bin)
    
    ix = [];
    if it == length(x_bin)
        ix = find(M(:,1)  >= x_bin(it));
    else
        ix = find(M(:,1)  >= x_bin(it) & M(:,1) < x_bin(it+1));
    end;
    
    c = zeros(length(y_bin),1);
    for jt = 1:length(y_bin)
        
        ix2 = [];
        if jt == length(y_bin)
            ix2 = find(M(ix,2)  >= y_bin(jt));
        else
            ix2 = find(M(ix,2)  >= y_bin(jt) & M(ix,2) < y_bin(jt+1));
        end;
        
        c(jt) = sum([length(ix2) 0]);%./length(y_bin);
        
        
    end;
    
    %smc = conv([zeros(1,length(c)) c' zeros(1,length(c))],gausswin(5),'same');
    smc = c;%smc(length(c)+1:end-length(c));
    
    x = smc./size(M,1);
    x(isnan(x)) = 0;
    N(:,it) = x;
    
end;
%%
hold on;
pcolor(x_bin,y_bin,log(N));
%contour(x_bin,y_bin,log(N),'w');
colormap hot;
shading interp;lighting phong;
axis tight;
axis xy;
cb = colorbar;
z = get(cb,'YLabel');
set(z,'String','PDF');
l = get(cb,'YTick');
l = [min(l) max(l)];
caxis(l);
set(cb,'YTick',l);
set(cb,'YTickLabel',round(exp(l).*10)./10);


return;