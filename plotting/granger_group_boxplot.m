clear gc; close all
for i = 1:7
    for side=1:2
a = gc_dist_sub_save{i,side};
gc_for(side,i,1) = mean(a{1});
gc_for(side,i,2) = mean(a{2});
a = gc_cd_dist_sub_save{i,side};
gc_bac(side,i,1) = mean(a{1});
gc_bac(side,i,2) = mean(a{2});
    end
end

ON = gc_for(:,:,1);
OFF = gc_for(:,:,2);
subplot(2,1,1);
barplot_N2(ON(:),OFF(:),'GC Backward','Granger'); %ylim([-1 0])

ON = gc_bac(:,:,1);
OFF = gc_bac(:,:,2);
subplot(2,1,2); 
barplot_N2(ON(:),OFF(:),'GC Forward','Granger'); %ylim([0 1.2])