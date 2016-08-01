
%%
p2d = '/bcbl/home/home_a-f/froux/chronset/thresholds/';

d1 = load([p2d,'greedy_optim_NP_data_BCN_28-Sep-2015.mat']);
d2 = load([p2d,'greedy_optim_SayWhen_18-Sep-2015.mat']);

d3 = load([p2d,'greedy_optim_NP_data_BCN_18-Jun-2016.mat']);
d4 = load([p2d,'greedy_optim_NP_data_SayWhen_18-Jun-2016.mat']);
%%
m1 = zeros(4,1);
for it = 1:4
    
    eval(['d',num2str(it),'.optim_data.hist_e(d',num2str(it),'.optim_data.hist_e==0)=NaN;']);
    
    ind = eval(['find(d',num2str(it),'.optim_data.hist_e == min(min(d',num2str(it),'.optim_data.hist_e)))']);
    [i1,i2] = eval(['ind2sub(size(d',num2str(it),'.optim_data.hist_e),ind)']);
    
    i1 = min(unique(i1));
    i2 = min(unique(i2));
    
    m1(it) = eval(['d',num2str(it),'.optim_data.hist_e(i1,i2)']);
    
end;

%%
m2 = zeros(4,1);
for it = 1:4
    
    eval(['d',num2str(it),'.optim_data.test_e(d',num2str(it),'.optim_data.test_e==0)=NaN;']);
    
    ind = eval(['find(d',num2str(it),'.optim_data.test_e == min(min(d',num2str(it),'.optim_data.test_e)))']);
    [i1,i2] = eval(['ind2sub(size(d',num2str(it),'.optim_data.test_e),ind)']);
    
    i1 = min(unique(i1));
    i2 = min(unique(i2));
    
    m2(it) = eval(['d',num2str(it),'.optim_data.test_e(i1,i2)']);
    
end;
%%
cmp = {'BCNold' m1(1) m2(1);...
    'SWold' m1(2) m2(2);...    
    'BCNnew' m1(3) m2(3);...    
    'SWnew' m1(4) m2(4);...
}
