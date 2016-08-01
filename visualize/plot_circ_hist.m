function plot_circ_hist(xr,yr)
%%
%%
d = max(xr)*2;
r = d/2;
th = 0:pi/360:2*pi;
x_u = r*cos(th)+0;
y_u = r*sin(th)+0;

x_u = normalize_data(x_u);
y_u = normalize_data(y_u);

xr = normalize_data(xr);
yr = normalize_data(yr);

return;