%%
clear;
clc;
%%
dataset = 'SW';%dataset to visualize
n = 1;%number of files to load
%%
switch dataset
    case 'SW'
        p2d = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/feature_data/';
        p2d2 = '/bcbl/home/home_a-f/froux/chronset/data/SayWhen/feature_data/Jun2016/';
        files = dir([p2d2,'*.mat']);
    case 'BCN'
        p2d = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/';
        p2d2 = '/bcbl/home/home_a-f/froux/chronset/data/BCN/feature_data/Jun2016/';
        files = dir([p2d2,'*.mat']);
end;
%%
if matlabpool('size') == 0
    matlabpool local;
end;
%%
rand_sel = randperm(length(files));
rand_sel = rand_sel(1:n);
files = files(rand_sel);

dat1 = cell(1,length(files));
dat2 = cell(1,length(files));
parfor it = 1:length(files)
    
    chck = dir([p2d,files(it).name(1:max(regexp(files(it).name,'_'))-1),'.mat']);
    
    dum = load([p2d,chck.name]); 
    dat1{it} = dum.savedata;
    
    dum = load([p2d2,files(it).name]);    
    dat2{it} = dum.savedata;
    
end;
%%
p1 = [.01:.003:.013];
p2 = [500:250:750];

figure;
for it = 1:length(dat2)
    k = 0;
    for kt = 1:length(p1)
        for lt = 1:length(p2)
            
            k = k+1;
            subplot(5,5,k);
            hold on;
            plot(dat1{it}.t,dat1{it}.features{1},'bo','LineWidth',3);
            plot(dat2{it}{kt,lt}.t,dat2{it}{kt,lt}.features{1},'r.','LineWidth',3);
            legend('before 06/16','after 06/16');
            set(gcf,'Color','w');
            title({dataset;['tw =',num2str(p1(kt)),'; sm =',num2str(p2(lt))]});
            
        end;
    end;
    %pause;
    %clf;
end;