function visualizeation_ds_comparions(ds)

[b1,bint,r,rint,stats1] = regress(ds.g,[ones(length(ds.g),1) ds.lds]);
[b2,bint,r,rint,stats2] = regress(ds.g,[ones(length(ds.g),1) ds.ds1]);
[b3,bint,r,rint,stats3] = regress(ds.g,[ones(length(ds.g),1) ds.ds2]);
b1
b2
b3
yp = b1(1)+b1(2)*ds.lds;

subplot(221);
hold on;
plot(ds.lds,yp,'r','LineSmoothing','on','LineWidth',3);
plot(ds.lds,ds.g,'o','Color',[.75 .75 .75],'MarkerFaceColor','k');
axis tight;


stats1(1)
stats2(1)
stats3(1)

subplot(223);

a = abs(ds.lds-ds.g);

b = linspace(0,200,500);


a1 = abs(ds.lds-ds.g);
a2 = abs(ds.ds1-ds.g);
a3 = abs(ds.ds2-ds.g);

[mean(a1) mean(a2) mean(a3)]

[n1,x1] = hist(a1,b);
[n2,x2] = hist(a2,b);
[n3,x3] = hist(a3,b);

hold on;
plot(x3,n3,'r-','LineSmoothing','on','LineWidth',3);%SW
plot(x2,n2,'b-','LineSmoothing','on','LineWidth',3);%CH
plot(x1,n1,'k-','LineSmoothing','on','LineWidth',3);%MM

axis tight;

subplot(224);

a = (ds.lds-ds.g);
b = linspace(-200,200,500);

a4 = (ds.lds-ds.g);
a5 = (ds.ds1-ds.g);
a6 = (ds.ds2-ds.g);

[mean(a4) mean(a5) mean(a6)]

[n1,x1] = hist(a4,b);
[n2,x2] = hist(a5,b);
[n3,x3] = hist(a6,b);



[n1,x1] = hist((ds.lds-ds.g),b);
[n2,x2] = hist((ds.ds1-ds.g),b);
[n3,x3] = hist((ds.ds2-ds.g),b);

hold on;

plot(x3,cumsum(n3)./sum(n3),'r-','LineSmoothing','on','LineWidth',3);
plot(x2,cumsum(n2)./sum(n2),'b-','LineSmoothing','on','LineWidth',3);
plot(x1,cumsum(n1)./sum(n1),'k-','LineSmoothing','on','LineWidth',3);

axis tight;

return
