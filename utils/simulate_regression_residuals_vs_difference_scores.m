%%
x = 500+randn(15,1);

X1 = x+10;
Y1 = x; 

x = 500+randn(15,1);
X2 = x;
Y2 = x+10;

x = 500+randn(15,1);
X3 = x;
Y3 = x+(10+1.*randn(15,1));

figure;
subplot(121);
a = gca;
hold on;
[b1,bint,r1,rint,stats] = regress(Y1,[ones(size(X1)) X1]);
yp1 = b1(1)+b1(2)*X1;
plot(X1,yp1,'r');
plot(X1,Y1,'bo','MarkerFaceColor','b');

subplot(122);
a = [a gca];
hold on;

[b2,bint,r2,rint,stats] = regress(Y2,[ones(size(X2)) X2]);
yp2 = b2(1)+b2(2)*X2;
plot(X2,yp2,'r');
plot(X2,Y2,'bo','MarkerFaceColor','b');

d1 = X1-Y1;
d2 = X2-Y2;
d3 = X3-Y3;

abs(mean(d1))
mean(abs(r1))


mean(abs(d2))
mean(abs(r2))


mean(abs(d3))
mean(abs(r3))

for it = 1:length(a)
    xlabel(a(it),'Automatic [ms]');
    ylabel(a(it),'Manual [ms]');
end;
set(a,'Box','off');
set(gcf,'Color','w');