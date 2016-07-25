%%
bw = .65;

figure;
subplot(211);
a = gca;
X = [hist_dat.onBCN1(:,2) hist_dat.onBCN2(:,2) hist_dat.onBCN3(:,2) hist_dat.onBCN1(:,1) hist_dat.onBCN2(:,1) hist_dat.onBCN3(:,1)];

X(:,4) = X(:,4)./sum(X(:,1));
X(:,5) = X(:,5)./sum(X(:,2));
X(:,6) = X(:,6)./sum(X(:,3));

X(:,1) = X(:,1)./sum(X(:,1));
X(:,2) = X(:,2)./sum(X(:,2));
X(:,3) = X(:,3)./sum(X(:,3));

hold on;

x = 1:3:size(X,1)*3;
k = 0;
for it = 1:size(X,1)    
    bar(x(it)+k,X(it,1),bw);
    k = k+1;
end;
h = [];
h = findobj('Type','patch');
set(h,'FaceColor','r');

x = 2:3:size(X,1)*3;
k = 0;
for it = 1:size(X,1)    
    bar(x(it)+k,X(it,2),bw);
    k = k+1;
end;
h = [];
h = findobj('Type','patch');
h = h(1:length(x));
set(h,'FaceColor',[.5 .5 .5]);

x = 3:3:size(X,1)*3;
k = 0;
for it = 1:size(X,1)    
    bar(x(it)+k,X(it,3),bw);
    k = k+1;
end;
h = [];
h = findobj('Type','patch');
h = h(1:length(x));
set(h,'FaceColor',[0 0 0]);

% k = 0;
% c = 0;
% for it = 1:3:size(X,1)*3
%     c = c+1;
%     plot([it it]+k,[X(c,1) X(c,1)+X(c,4)],'k');
%     plot([it-.3 it+.3]+k,[X(c,1)+X(c,4) X(c,1)+X(c,4)],'k');
%     k = k+1;
% end;
% 
% k = 0;
% c = 0;
% for it = 2:3:size(X,1)*3
%     c = c+1;
%     plot([it it]+k,[X(c,2) X(c,2)+X(c,5)],'k');
%     plot([it-.3 it+.3]+k,[X(c,2)+X(c,5) X(c,2)+X(c,5)],'k');
%     k = k+1;
% end;
% 
% k = 0;
% c = 0;
% for it = 3:3:size(X,1)*3
%     c = c+1;
%     plot([it it]+k,[X(c,3) X(c,3)+X(c,6)],'k');
%     plot([it-.3 it+.3]+k,[X(c,3)+X(c,6) X(c,3)+X(c,6)],'k');
%     k = k+1;
% end;

x = 2:4:size(X,1)*4;
x(2:end) = x(2:end)+1;
bs = [5:5:25 50:25:150];
bs = sort([-bs bs]);

set(gca,'XTick',x);
set(gca,'XTickLabel',{'Less' bs 'More'});

axis tight;

subplot(212);
a = [a gca];
X = [hist_dat.onSW1(:,2) hist_dat.onSW2(:,2) hist_dat.onSW3(:,2) hist_dat.onSW4(:,2) hist_dat.onSW1(:,1) hist_dat.onSW2(:,1) hist_dat.onSW3(:,1) hist_dat.onSW4(:,1)];

bw = .9;

X(:,1) = X(:,1)./sum(X(:,1));
X(:,2) = X(:,2)./sum(X(:,2));
X(:,3) = X(:,3)./sum(X(:,3));
X(:,4) = X(:,4)./sum(X(:,4));


hold on;

x = 1:4:size(X,1)*4;
k = 0;
for it = 1:size(X,1)
    bar(x(it)+k,X(it,1),bw);
    k = k+1;    
end;
h = [];
h = findobj('Type','patch');
h = h(1:length(x));
set(h,'FaceColor','r');

x = 2:4:size(X,1)*4;
k = 0;
for it = 1:size(X,1)
    bar(x(it)+k,X(it,2),bw);
    k = k+1;    
end;
h = [];
h = findobj('Type','patch');
h = h(1:length(x));
set(h,'FaceColor',[.5 .5 .5]);

x = 3:4:size(X,1)*4;
k = 0;
for it = 1:size(X,1)
    bar(x(it)+k,X(it,3),bw);
    k = k+1;    
end;
h = [];
h = findobj('Type','patch');
h = h(1:length(x));
set(h,'FaceColor',[0 0 0]);

x = 4:4:size(X,1)*4;
k = 0;
for it = 1:size(X,1)
    bar(x(it)+k,X(it,4),bw);
    k = k+1;    
end;
h = [];
h = findobj('Type','patch');
h = h(1:length(x));
set(h,'FaceColor',[1 1 1]);

% k = 0;
% c = 0;
% for it = 1:4:size(X,1)*4
%     c = c+1;
%     plot([it it]+k,[X(c,1) X(c,1)+X(c,5)],'k');
%     plot([it-.3 it+.3]+k,[X(c,1)+X(c,5) X(c,1)+X(c,5)],'k');
%     k = k+1;
% end;
% 
% k = 0;
% c = 0;
% for it = 2:4:size(X,1)*4
%     c = c+1;
%     plot([it it]+k,[X(c,2) X(c,2)+X(c,6)],'k');
%     plot([it-.3 it+.3]+k,[X(c,2)+X(c,6) X(c,2)+X(c,6)],'k');
%     k = k+1;
% end;
% 
% k = 0;
% c = 0;
% for it = 3:4:size(X,1)*4
%     c = c+1;
%     plot([it it]+k,[X(c,3) X(c,3)+X(c,7)],'k');
%     plot([it-.3 it+.3]+k,[X(c,3)+X(c,7) X(c,3)+X(c,7)],'k');
%     k = k+1;
% end;
% 
% k = 0;
% c = 0;
% for it = 4:4:size(X,1)*4
%     c = c+1;
%     plot([it it]+k,[X(c,4) X(c,4)+X(c,8)],'k');
%     plot([it-.3 it+.3]+k,[X(c,4)+X(c,8) X(c,4)+X(c,8)],'k');
%     k = k+1;
% end;

axis tight;

x = 2:4:size(X,1)*4;

c = 0;
for it = 1:length(x)
    x(it) = x(it)+c;
    c = c+1;
end;
bs = [5:5:25 50:25:150];
bs = sort([-bs bs]);

set(gca,'XTick',x);
set(gca,'XTickLabel',{'Less' bs 'More'});


for it = 1:length(a)
    xlabel(a(it),'Absolute dfierence automatic-manuals [ms]');
    ylabel(a(it),'Proportion of trials [%]');
end;

set(gcf,'Color','w');