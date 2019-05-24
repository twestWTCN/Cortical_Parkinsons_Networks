function [] = barplot_coh_groups(datapathr,highbetacoh)
% % subplot(1,2,1)
% % sideCoh(:,:,2,1) = vertcat(highbetacoh{2,1,:}); % Left ON
% % sideCoh(:,:,1,1) = vertcat(highbetacoh{1,1,:}); % Left OFF
% % sideCoh(:,:,2,2) = vertcat(highbetacoh{2,2,:}); % Right ON
% % sideCoh(:,:,1,2) = vertcat(highbetacoh{1,2,:}); % Right OFF
% % 
% % boxplot([squeeze(sideCoh(:,1,1,1)) squeeze(sideCoh(:,1,2,1)) squeeze(sideCoh(:,1,1,2)) squeeze(sideCoh(:,1,2,2))],'labels',{'Left OFF','Left ON','Right OFF','Right ON'},...
% %     'BoxStyle','filled','Widths',0.8)
% % %%
% % a = get(get(gca,'children'),'children');   % Get the handles of all the objects
% % t = get(a,'tag');   % List the names of all the objects 
% % idx=strcmpi(t,'box');  % Find Box objects
% % boxes=a(idx);          % Get the children you need
% % set(boxes,'linewidth',35); % Set width
% % idx=strcmpi(t,'Whisker');  % Find Box objects
% % whisker=a(idx);          % Get the children you need
% % set(whisker,'LineWidth',2); % Set width
% % idx=strcmpi(t,'Median');  % Find Box objects
% % whisker=a(idx);          % Get the children you need
% % set(whisker,'LineWidth',2); % Set width
% % %%
% % grid on; ylabel('Coherence'); title('Hemi Analyses');  
% % % tests
% % % Left ON vs OFF
% % [h pv(1)] = ttest2(squeeze(sideCoh(:,1,1,1)),squeeze(sideCoh(:,1,2,1)));
% % % Right ON vs OFF
% % [h pv(2)] = ttest2(squeeze(sideCoh(:,1,1,2)),squeeze(sideCoh(:,1,2,2)));
% % % Left ON vs Right ON
% % [h pv(3)] = ttest2(squeeze(sideCoh(:,1,1,1)),squeeze(sideCoh(:,1,1,2)));
% % % Left OFF vs Right OFF
% % [h pv(4)] = ttest2(squeeze(sideCoh(:,1,2,1)),squeeze(sideCoh(:,1,2,2))); 
% % H=sigstar({{'Left OFF','Left ON'},{'Right OFF','Right ON'},{'Left OFF','Right ON'},{'Left OFF','Right ON'}},pv); 
% % ylim([0 0.7])

pv = [];
% subplot(1,2,2)
Coh(:,:,1) = vertcat(highbetacoh{1,:,:});
Coh(:,:,2) = vertcat(highbetacoh{2,:,:});
[h pv] = ttest2(squeeze(Coh(:,1,1)),squeeze(Coh(:,1,2)));
boxplot([Coh(:,1,1) Coh(:,1,2)],'labels',{'Total OFF','Total ON'},...
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
grid on; ylabel('Coherence'); title('Bilateral Analyses') 
H=sigstar({{'Total OFF','Total ON'}},pv); ylim([0 0.5])
set(gcf,'Position',[1217         479         400         350])

