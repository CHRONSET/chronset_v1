function [T1,T2] = monte_carlo_simulation_samples_size_vs_tstat3(M,SD,AN,params)
%%
clear;
clc;
%%
if matlabpool('size')==0
    matlabpool 128;
end;
%% parameters for simulation
if nargin == 0
    M = 15:15:30;% mean of sample
    SD = 50:25:125;% sd of sample
    AN = [90];%chronset uncertainty
    params.N = 5:50;% number of participants (sample size)
    params.n = 100;% number of obs per participants
    params.n_iter = 20e4;%number of sample draws
end;
%% preallocate matrix with t-values 
T1 = zeros(length(M),length(SD),length(params.N));
T2 = zeros(length(M),length(SD),length(params.N),length(AN));
%% do the simulation
for mt = 1:length(M) %loop over sample means
    for nt = 1:length(SD)% loop over sample SD
                        
        for jt = 1:length(params.N)% loop over samples sizes
        
            t_val1 = zeros(params.n_iter,1);
            t_val2 = zeros(params.n_iter,length(AN));
            
            tic;
            parfor it = 1:params.n_iter%loop over number of sample draws
                
                %pre-allocate
                x1 = zeros(params.N(jt),params.n);
                x2 = zeros(params.N(jt),params.n);
                
                x3 = zeros(params.N(jt),length(AN),params.n);
                x4 = zeros(params.N(jt),length(AN),params.n);
                
                for kt = 1:params.N(jt)% simulate n observations for each individual participant
                    %empirical
                    
                    %by adding the noise to EXACTLY the same vectors, we will
                    %better control for variability between the two distributions
                    v1 = 0+SD(nt).*randn(1,params.n);
                    v2 = M(mt)+SD(nt).*randn(1,params.n);

                    x1(kt,:) = v1;
                    x2(kt,:) = v2;
                    
                    for ot = 1:length(AN)% loop over additive noise
                        % add different levels of noise
                        x3(kt,ot,:) = v1+AN(ot).*randn(1,params.n);
                        x4(kt,ot,:) = v2+AN(ot).*randn(1,params.n);

                    end;     
                    
                end;
                
                % average over n observations for each participant
                x1 = mean(x1,2);
                x2 = mean(x2,2);
                x3 = squeeze(mean(x3,3));
                x4 = squeeze(mean(x4,3));
                % each participant is now an average over n-observations                                
                
                % dependent samples t-test of the mean against zero
                d = x2-x1;
                [~,~,~,stats] = ttest(d,0);% empirical parameters
                t_val1(it) = stats.tstat;
                
                d = x4-x3;
                [~,~,~,stats] = ttest(d,0);% noise added
                t_val2(it,:) = stats.tstat;
                               
            end;    
            toc;
            
            [trsh] = tinv(0.97250,params.N(jt)-1);
            
            n1 = length(find(t_val1 > trsh));
            T1(mt,nt,jt) = n1;
            
            n2 = zeros(size(t_val2,2),1);
            for ot = 1:size(t_val2,2)
                n2(ot) = length(find(t_val2(:,ot) > trsh));
            end;
            T2(mt,nt,jt,:) = n2';
            
        end;
    end;
end;
save(['/bcbl/home/home_a-f/froux/chronset/thresholds/mc_confidence_estimates_',date,'.mat'],'T1','T2');

matlabpool close;

return;

% %% plot results of simulation
% x1 = zeros(N(end),n);
% x2 = zeros(N(end),n);
% x3 = zeros(N(end),n);
% x4 = zeros(N(end),n);
% for it = 1:N(end)
%     x1(it,:) = M1+SD1*randn(1,n);
%     x2(it,:) = M2+SD2*randn(1,n);
%     x3(it,:) = M1+SD2*randn(1,n)+AN*randn(1,n);
%     x4(it,:) = M1+SD2*randn(1,n)+AN*randn(1,n);
% end;
%
% y1 = mean(x1,2);
% y2 = mean(x2,2);
% y3 = mean(x3,2);
% y4 = mean(x4,2);
%
% [n1,x1] = hist(y1,min(y1)-2:.5:max(y1)+2);
% [n2,x2] = hist(y2,min(y2)-2:.5:max(y2)+2);
% [n3,x3] = hist(y3,min(y3)-2:.5:max(y3)+2);
% [n4,x4] = hist(y4,min(y4)-2:.5:max(y4)+2);
%
% figure;
% subplot(221);
% a = gca;
% hold on;
% bar(x1,n1,.95);
% bar(x2,n2,.95,'FaceColor','r');
%
% subplot(222);
% a = [a gca];
% hold on;
% bar(x3,n3,.95);
% bar(x4,n4,.95,'FaceColor','r');
%
% axis(a,'tight');
% set(a(1),'YLim',[0 max([n1 n2])]);
% set(a(2),'YLim',[0 max([n3 n4])]);
%
% set(a(1),'YTick',[0 max([n1 n2])]);
% set(a(2),'YTick',[0 max([n3 n4])]);
%
%
% z = norminv(1-.1,0,1); p = 1-normpdf(z,0,1);[z p];
%
% x = N(1:5:end);
% y1 = mean(t1(:,1:5:end),1);
% y2 = z*std(t1(:,1:5:end),1);
%
% y3 = mean(t2(:,1:5:end),1);
% y4 = z*std(t2(:,1:5:end),1);
%
% %DO WE NOT ACTUALLY WANT SOMETHING LIKE THE 20th percentile for the
% %bottom?
% d1 =min(find(sign((y1-y2)-z)==1));
% d2 =min(find(sign((y3-y4)-z)==1));
%
%
%
% %%Final plot of t-values for both tests
% subplot(2,2,3);
% a = [a gca];
% hold on;
% plot([N(1) N(end)],[2 2],'r--','LineWidth',3);
%
%
% %plot([n(d1) n(d1)],[0 0],'r^','LineWidth',3);
% %plot([n(d2) n(d2)],[0 0],'b^','LineWidth',3);
%
% plot(x,y1,'-','MarkerFaceColor',[.75 .75 .75],'LineWidth',3);
% for jt = 1:length(x)
%
%     %plot([x(jt) x(jt)],[y1(jt)-y2(jt) y1(jt)],'k','LineWidth',3);
%     %plot([x(jt)-3.5 x(jt)+3.5],[y1(jt)-y2(jt) y1(jt)-y2(jt)],'k','LineWidth',3);
% end;
%
% plot(x,y3,'-.','MarkerFaceColor',[.75 .75 .75],'LineWidth',3);
% for jt = 1:length(x)
%
%     %plot([x(jt) x(jt)],[y3(jt)-y4(jt) y3(jt)],'k','LineWidth',3);
%     %plot([x(jt)-3.5 x(jt)+3.5],[y3(jt)-y4(jt) y3(jt)-y4(jt)],'k','LineWidth',3);
% end;
%
% axis(a,'tight');
% set(a,'Fontsiz',14);
% xlabel(a(end),'Sample size','Fontsize',14);
% ylabel(a(end),'t-statistic [a.u.]','Fontsize',14);
%
% %set(gca,'YTick',[5:5:max(max((t)))]);
% %set(gca,'XTick',[100:100:1000]);
%
%
% set(gcf,'Color','w');
%
% %print out the maximum difference in the vectors for the figure caption
% disp('Max difference between y1 (empirical t) and y3 (chron_t)');
% disp(max(y1-y3))

