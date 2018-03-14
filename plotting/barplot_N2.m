function [] = barplot_N2(A,B,ylab,tit)
boxplot([A B],'labels',{'ON','OFF'},...
    'BoxStyle','filled','Widths',0.8);
h=findobj(gca,'tag','Outliers'); 
    delete(h)
%%
a = get(get(gca,'children'),'children');   % Get the handles of all the objects
t = get(a,'tag');   % List the names of all the objects 
idx=strcmpi(t,'box');  % Find Box objects
boxes=a(idx);          % Get the children you need
set(boxes,'linewidth',35); % Set width
idx=strcmpi(t,'Whisker');  % Find Box objects
whisker=a(idx);          % Get the children you need
set(whisker,'LineWidth',2); % Set width
idx=strcmpi(t,'Median');  % Find Box objects
whisker=a(idx);          % Get the children you need
set(whisker,'LineWidth',2); % Set width
%%
grid on; ylabel(ylab); title(tit);  
% tests
% ON vs OFF
[h pv(1)] = ttest2(A,B);
% [h pv(1)] = vartest2(A,B)
H=sigstar({{'ON','OFF'}},pv);
ylim(get(gca,'YLim').*[1 1.05])

