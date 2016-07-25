%%
[u,s,v] = svd(Y);
recon = u(:,1)*s(1,1)*v(:,1)';
D1 = diff(Y,[],2);
D2 = diff(recon,[],2);
recon = mean(recon,2);


[b1,bint1,r1,rint1,stats1] = regress(recon,[ones(size(recon)) Y(:,1)]);
[b2,bint2,r2,rint2,stats2] = regress(recon,[ones(size(recon)) on1]);
[b3,bint3,r3,rint3,stats3] = regress(recon,[ones(size(recon)) on2]);
[b4,bint4,r4,rint4,stats4] = regress(recon,[ones(size(recon)) on3]);


yp1 = b1(1)+b1(2).*Y(:,1);
yp2 = b2(1)+b2(2).*on1;
yp3 = b3(1)+b3(2).*on2;
yp4 = b4(1)+b4(2).*on3;



R1 = recon-yp1;
R2 = recon-yp2;
R3 = recon-yp3;
R4 = recon-yp4;

R = [R1 R2 R3 R4];
m1 = min(min(R));
m2 = max(max(R));

[n1,x1] = hist(R1,linspace(m1,m2,250));
[n2,x2] = hist(R2,linspace(m1,m2,250));
[n3,x3] = hist(R3,linspace(m1,m2,250));
[n4,x4] = hist(R4,linspace(m1,m2,250));

g1 = n1./sum(n1);
g2 = n2./sum(n2);
g3 = n3./sum(n3);
g4 = n4./sum(n4);

figure;

a = zeros(8,1);
for it = 1:8
    subplot(4,2,it);
    a(it) = gca;
    hold(a(it),'on');
end;

subplot(4,2,1);plot(sort(D1),'b-');plot(sort(D2),'r-');


subplot(4,2,2);
plot(Y(:,1),yp1,'Color',[.25 .25 .75]);
plot(on1,yp2,'Color',[.25 .5 .5]);
plot(on2,yp3,'Color',[.5 0 .75]);
plot(on3,yp4,'Color',[.75 .75 .75]);

%plot(Y(:,1),recon,'k.');


subplot(4,2,3);plot(yp1,zeros(length(r),1),'r-');plot(yp1,r1,'k.');title(num2str(sqrt(mean(r1.^2))));

subplot(4,2,4);plot(yp2,zeros(length(r),1),'r-');plot(yp2,r2,'k.');title(num2str(sqrt(mean(r2.^2))));

subplot(4,2,5);plot(yp3,zeros(length(r),1),'r-');plot(yp3,r3,'k.');title(num2str(sqrt(mean(r3.^2))));

subplot(4,2,6);plot(yp4,zeros(length(r),1),'r-');plot(yp4,r4,'k.');title(num2str(sqrt(mean(r4.^2))));


subplot(4,2,7);
plot(x1,cumsum(g1),'Color',[.25 .25 .75]);
plot(x1,cumsum(g2),'Color',[.25 .5 .5]);
plot(x1,cumsum(g3),'Color',[.5 0 .75]);
plot(x1,cumsum(g4),'Color',[.75 .75 .75]);

axis(a,'tight');
set(a(3:6),'YLim',[m1 m2]);
set(a(7),'YLim',[0 1.05]);

xlabel(a(2),'Measured onset [ms]');
ylabel(a(2),'\Lambda [ms]');

for it = 3:length(a)
    
    xlabel(a(it),'Predicted onset [ms]');
    ylabel(a(it),'Residual [ms]');
end;



set(gcf,'Color','w');
%%
nLogLGradFun = @(theta) deal(-sum(-gammaln(theta(1)) - ...
    theta(1)*log(theta(2) + x) + (theta(1)-1)*log(y) -...
    y./(theta(2)+x)),...
    -[sum(-psi(theta(1))+log(y./(theta(2)+x)));...
    sum(1./(theta(2)+x).*(y./(theta(2)+x)-theta(1)))]);

x = Y(:,1);
y = recon;

theta0 = randn(2,1);
uLB = [0 -min(x)];
uUB = [Inf Inf];
options = optimset('Algorithm','active-set',...
    'TolFun',1e-10,'Display','off','GradObj','on');

[uMLE,uLogL] = fmincon(nLogLGradFun,theta0,[],[],[],[],uLB,uUB,[],options);
uLogL = -uLogL;