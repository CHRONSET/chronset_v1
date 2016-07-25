%%
clear;
clc;
restoredefaultpath;
%% parameters for simulation

M1 = 0;% these are the empirical means and SDs
SD1 = 11;

M2 = 10;
SD2 = 11;

M3 = 0;%these are the noise added means and SDs
SD3 = 30;
M4 = 10;
SD4 = 30;

N = 10:400;% number of participants (sample size)
n = 100;% number of observations per participants
n_iter = 100;% number of simulations per sample size
%% preallocate t-values and confidence intervals
t1 = zeros(n_iter,length(N));
t2 = zeros(n_iter,length(N));
c1 = zeros(n_iter,length(N));
c2 = zeros(n_iter,length(N));
%% do the simulation
for it = 1:n_iter%loop over number of emulations
    tic;
    
    %pre-allocate
    t_val1 = zeros(length(N),1);
    t_val2 = zeros(length(N),1);
    ci1 = zeros(length(N),1);
    ci2 = zeros(length(N),1);
    
    for jt = 1:length(N)% loop over samples sizes        
        
        %pre-allocate
        x1 = zeros(N(jt),n);
        x2 = zeros(N(jt),n);
        x3 = zeros(N(jt),n);
        x4 = zeros(N(jt),n);
        
        for kt = 1:N(jt)% simulate n observations for each individual participant
            %empirical
            x1(kt,:) = M1+SD1*randn(1,n);
            x2(kt,:) = M2+SD2*randn(1,n);
            % noise added
            x3(kt,:) = M3+SD3*randn(1,n);            
            x4(kt,:) = M4+SD4*randn(1,n);
            
        end;
        
        % average over n observations for each participant
        x1 = mean(x1,2);
        x2 = mean(x2,2);
        x3 = mean(x3,2);
        x4 = mean(x4,2);
        % each participant is now an average over n-observations
        
        % independent samples t-test of the means between conditions
        [h,p,ci,stats] = ttest2(x2',x1');% empirical parameters
        t_val1(jt) = stats.tstat;
        ci1(jt) = min(ci);
        
        [h,p,ci,stats] = ttest2(x4',x3');% noise added parameterts
        t_val2(jt) = stats.tstat;
        ci2(jt) = min(ci);
        
        
    end;
    
    t1(it,:) = t_val1;
    t2(it,:) = t_val2;
    c1(it,:) = ci1;
    c2(it,:) = ci2;
        
    toc;
end;
%% plot results of simulation
x1 = zeros(N(end),n);
x2 = zeros(N(end),n);
x3 = zeros(N(end),n);
x4 = zeros(N(end),n);
for it = 1:N(end)
    x1(it,:) = M1+SD1*randn(1,n);
    x2(it,:) = M2+SD2*randn(1,n);
    x3(it,:) = M3+SD3*randn(1,n);
    x4(it,:) = M4+SD4*randn(1,n);
end;

y1 = mean(x1,2);
y2 = mean(x2,2);
y3 = mean(x3,2);
y4 = mean(x4,2);

[n1,x1] = hist(y1,min(y1)-2:.5:max(y1)+2);
[n2,x2] = hist(y2,min(y2)-2:.5:max(y2)+2);
[n3,x3] = hist(y3,min(y3)-2:.5:max(y3)+2);
[n4,x4] = hist(y4,min(y4)-2:.5:max(y4)+2);

figure;
subplot(221);
a = gca;
hold on;
bar(x1,n1,.95);
bar(x2,n2,.95,'FaceColor','r');

subplot(222);
a = [a gca];
hold on;
bar(x3,n3,.95);
bar(x4,n4,.95,'FaceColor','r');

axis(a,'tight');
set(a(1),'YLim',[0 max([n1 n2])]);
set(a(2),'YLim',[0 max([n3 n4])]);

set(a(1),'YTick',[0 max([n1 n2])]);
set(a(2),'YTick',[0 max([n3 n4])]);


z = norminv(1-.1,0,1); p = 1-normpdf(z,0,1);[z p];

x = N(1:5:end);
y1 = mean(t1(:,1:5:end),1);
y2 = z*std(t1(:,1:5:end),1);

y3 = mean(t2(:,1:5:end),1);
y4 = z*std(t2(:,1:5:end),1);

d1 =min(find(sign((y1-y2)-z)==1));
d2 =min(find(sign((y3-y4)-z)==1));

subplot(2,2,3);
a = [a gca];
hold on;
plot([N(1) N(end)],[2 2],'r--','LineWidth',3);


%plot([n(d1) n(d1)],[0 0],'r^','LineWidth',3);
%plot([n(d2) n(d2)],[0 0],'b^','LineWidth',3);

plot(x,y1,'ko','MarkerFaceColor',[.75 .75 .75],'LineWidth',3);
for jt = 1:length(x)
    
    plot([x(jt) x(jt)],[y1(jt)-y2(jt) y1(jt)],'k','LineWidth',3);
    plot([x(jt)-3.5 x(jt)+3.5],[y1(jt)-y2(jt) y1(jt)-y2(jt)],'k','LineWidth',3);
end;

plot(x,y3,'ko','MarkerFaceColor',[.75 .75 .75],'LineWidth',3);
for jt = 1:length(x)
    
    plot([x(jt) x(jt)],[y3(jt)-y4(jt) y3(jt)],'k','LineWidth',3);
    plot([x(jt)-3.5 x(jt)+3.5],[y3(jt)-y4(jt) y3(jt)-y4(jt)],'k','LineWidth',3);
end;

axis(a,'tight');
set(a,'Fontsiz',14);
xlabel(a(end),'Sample size','Fontsize',14);
ylabel(a(end),'t-statistic [a.u.]','Fontsize',14);

%set(gca,'YTick',[5:5:max(max((t)))]);
%set(gca,'XTick',[100:100:1000]);


set(gcf,'Color','w');