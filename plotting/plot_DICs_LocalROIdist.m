function plot_DICs_LocalROIdist(R)
%%%
% This function will plot the distribution of the local DICs peaks found
% within the global prespecified ROIs from Litvak
%%%
close all; clear
R = makeHeader_SubCort_Cort_Networks();
scatleg = {'bo','bx';'ro','rx'};
for band = [1 3] %:numel(R.bandname)
    % Pial Surface in MNI
    %     figure; set(gcf,'color','w');
    %     load([R.datapathr 'template_MRI\surface_pial_both.mat'])
    %         ft_plot_mesh(mesh);
    %         hold on
    %         a = gca;
    %         a.Children(1).FaceAlpha = 0.3;
    % OR grid in dim space
    figure
    load([R.datapathr 'template_MRI/standard_sourcemodel3d10mm']);
    template_grid = sourcemodel;
    clear sourcemodel;
    mask =  reshape(template_grid.inside,template_grid.dim)>0;
    hold on
    fv = isosurface(mask,0);
    patch(fv,'FaceColor',[.0 .0 1],'EdgeColor',[0 0 0],'FaceAlpha',0.05,'LineWidth',0.06);
    view([30 30])
    
    for sub = 1:numel(R.subname)
        load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_DICSv2_source_channelselectioninfo_' R.bandname{band}])
        for cond = 1:2
            [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
            [~,~,nrepOFF,~] = data_fileguide(R.subname{sub},1);
            for nr = 1:nrep
                if nrepOFF<nrep
                    nroff = 1;
                else
                    nroff = nr;
                end
                
                for side = 1:2
                    switch R.ipsicon
                        case 'ipsi'
                            seldet = STNselection{side,1,nrepOFF,cond}; % selected for ipsi % Use OFF
                        case 'contra'
                            seldet = STNselection{side,2,nrepOFF,cond}; % selected for contra % Use OFF
                    end
                    ROI = seldet{7};
                    ROI = [ROI(2) ROI(1) ROI(3)];
                    scatter3(ROI(1),ROI(2),ROI(3),125,scatleg{cond,side},'LineWidth',3);
                    hold on
                    disp([side nr cond sub])
                    view([-90 90]); axis equal
                end
            end
        end
    end
end
a = 1;
