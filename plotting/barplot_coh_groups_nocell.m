function [] = barplot_coh_groups_nocell(datapathr,cohmax)
subplot(1,2,1)
block = []; grouping = [];
x = squeeze(cohmax(1,:,:,1)); x(x==0) = []; %sideCoh(:,1,1,1) = x(:); % Left ON
block = [block x]; grouping = [grouping repmat(1,1,size(x,2))];
x = squeeze(cohmax(1,:,:,2));  x(x==0) = []; %sideCoh(:,1,2,1) = x(:); % Left OFF
block = [block x]; grouping = [grouping repmat(2,1,size(x,2))];
x = squeeze(cohmax(2,:,:,1));  x(x==0) = [];%sideCoh(:,1,1,2) = x(:); % Right ON
block = [block x]; grouping = [grouping repmat(3,1,size(x,2))];
x = squeeze(cohmax(2,:,:,2));  x(x==0) = [];%sideCoh(:,1,2,2) = x(:); % Right OFF
block = [block x]; grouping = [grouping repmat(4,1,size(x,2))];

boxplot(block,grouping,'labels',{'Left ON','Left OFF','Right ON','Right OFF'},...
    'BoxStyle','filled','Widths',0.8)
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
grid on; ylabel('Maximum WPLI'); title('Hemi Analyses');  
% tests
% Left ON vs OFF
[h pv(1)] = ttest2(block(grouping==1),block(grouping==2));
% Right ON vsbl OFF
[h pv(2)] = ttest2(block(grouping==3),block(grouping==4));
% Left ON vs Right ON
[h pv(3)] = ttest2(block(grouping==1),block(grouping==3));
% Left OFF vs Right OFF
[h pv(4)] = ttest2(block(grouping==2),block(grouping==4));
H=sigstar({{'Left ON','Left OFF'},{'Right ON','Right OFF'},{'Right ON','Left ON'},{'Right OFF','Left OFF'}},pv); ylim([0 0.3])

pv = []; block = []; grouping = []; Coh = {}
subplot(1,2,2)
x = squeeze(cohmax(:,:,:,1)); x(x==0) = []; 
block = [block x]; grouping = [grouping repmat(1,1,size(x(:),1))];
Coh{1} = x;
x = squeeze(cohmax(:,:,:,2)); x(x==0) = []; 
block = [block x]; grouping = [grouping repmat(2,1,size(x(:),1))];
Coh{2} = x;
[h pv] = ttest2(Coh{1},Coh{2});
boxplot(block,grouping,'labels',{'Total ON','Total OFF'},...
        'BoxStyle','filled','Widths',0.8)
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
grid on; ylabel('Maximum WPLI'); title('Bilateral Analyses') 
H=sigstar({{'Total ON','Total OFF'}},pv); ylim([0 0.3])
set(gcf,'Position',[512         540        1229         436])

