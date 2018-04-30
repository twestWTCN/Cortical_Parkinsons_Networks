function plot_phase_amp_analysis_PLIs_Cagnan(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
idd = '';
% Collate the points and then do paired ttest between bins? A lot of
% information is lost by cutting off in histograms. Either that or somehow
% increase the number of samples - perhaps include all of the STN channels?
% Maybe combine left and right into one histogram? Collapse all of one
% subjects sync segs into one histogram?? All should be sampling from the
% same distribution so build up better resolved image of underlying 2D PDF.
% Or look at doing a hotelling test between the collated data??

close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};

QX = 8 ; QY = 8; % 6 6
% QX = [-pi -pi/8 pi/8 pi]';
QX = linspace(-pi,pi,QX);
% load([R.datapathr 'subject_hbWPLI075'])
% subscreen = squeeze(sum(subject_hbcohscreen>R.PA.WPLIscreen)==2);
a = logspace(0.5,log10(150),4); logscalez = [-100 -50 -15 15 50 100]; %-40:10:40; %linspace(-50,50,QY); %logscalez = [-3.^(4:-1:2) 3.^(2:4)];
for band = 3%; %1:numel(R.bandname)
    for sub = 1:numel(R.subname)
        load([R.datapathr 'results\spectral\screen.mat'])
        for side =1:2
            if screen(side,sub) == 1
                
                for cond = 1:2
                    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
                    [datafileN,pp_mark,nrOFF,senscheck] = data_fileguide(R.subname{sub},1);
                    
                    clear relativePhi segL ampSeg HdistSeg tends
                    for nr = 1:nrep
                        if nr == 2;                            a = 1; end
                        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '.mat'])
                        relativePhi{nr} = wrapToPi(pA_pli_dist_save{cond,nr}-circ_mean(pA_pli_dist_save{2,nrOFF}') ); %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
                        segL{nr} =  segL_pli_dist_save{cond,nr}; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
                        ampSeg{nr} = amp_pli_dist_save{cond,nr}; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
                        HdistSeg{nr} = H_dist_save{cond,nr}(1,:);
                        tends(nr) = timevec{cond,nr}(end)-timevec{cond,nr}(1);
                    end
                    tendtot = sum(tends(:));
                    relativePhiCol = [relativePhi{:,:}];
                    segLCol = [segL{:,:}];
                    ampSegCol = [ampSeg{:,:}];
                    HdistSegCol = [HdistSeg{:,:}];
                    
                    phiBin = linspace(-pi,pi,10);
                    phiBinMid = phiBin(1:end-1)+((phiBin(2)-phiBin(1))/2);
                    ampBin = [];
                    for i = 1:length(phiBin)-1
                        binDat = ampSegCol(:,relativePhiCol>=phiBin(i) & relativePhiCol<=phiBin(i+1))';
                        ampBinMu(:,i) = nanmean(binDat);
                        ampBinSEM(:,i) = nanstd(binDat)/size(binDat,1);
                        
                        binSeg = segLCol(:,relativePhiCol>=phiBin(i) & relativePhiCol<=phiBin(i+1))';
                        segBinMu(:,i) = nanmean(binSeg);
                        segBinSEM(:,i) = nanstd(binSeg)/size(binSeg,1);
                    end
                    ampBinGroup{cond,sub} = ampBinMu;
                    segBinGroup{cond,sub} = segBinMu;
                    figure(1)
                    cmap = linspecer(3);
                    panlist = [1 3 5 ; 2 4 6]; ylimlist = {[-50 300];[-50 150];[-50 75]};
                    obs = {['PMC ' R.bandname{band}],['STN ' R.bandname{band}],['STN ' R.bandname{band-1}]};
                    for i = 1:3
                        subplot(3,2,panlist(cond,i))
                        [hl, hp] = boundedline(phiBinMid', ampBinMu(i,:)',ampBinSEM(i,:)','cmap',cmap(i,:));
                        if cond == 1; hl.LineStyle = '--'; end
                        hp.FaceAlpha = 0.4;
                        ylim(ylimlist{i})
                        title(['STN-PMC Phase vs ' obs{i}])
                        ylabel(['Amplification of ' obs{i}]); xlabel('Relative Phase')
                        grid on
                    end
                    
                    figure(2)
                    cmap = linspecer(5);
                    subplot(1,2,cond)
                    [hl, hp] = boundedline(phiBinMid', segBinMu(1,:)',segBinSEM(1,:)','cmap',cmap(5,:));
                    if cond == 1
                        hl.LineStyle = '--';
                    end
                    hp.FaceAlpha = 0.4;
                    ylim([0.15 0.75]); grid on
                    title(['STN-PMC Phase vs LB Frame Length'])
                    ylabel(['LB Segment Length (s)']); xlabel('Relative Phase')
                end
            end
        end
    end
    
    figure(1)
    set(gcf,'Position',[680          94         864        1004])
    for cond=1:2
        y =  horzcat(ampBinGroup{cond,:});
        for i = 1:3
            subplot(3,2,panlist(cond,i)); hold on
            a = plot(phiBinMid,nanmean(reshape(y(i,:),size(phiBinMid,2),[])'),'color',cmap(i,:),'LineWidth',2);
            if cond == 1; a.LineStyle = '--'; end
        end
    end
    annotation(gcf,'textbox',...
        [0.187 0.956 0.205 0.035],...
        'String',{'Control'},...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    annotation(gcf,'textbox',...
        [0.624 0.959 0.205 0.0349],...
        'String','Lesion',...
        'LineStyle','none',...
        'HorizontalAlignment','center',...
        'FontSize',18,...
        'FitBoxToText','off');
    
    figure(2)
    set(gcf,'Position',[681   761   862   336])
            cmap = linspecer(5);
    for cond=1:2
        y =  vertcat(segBinGroup{cond,:});
            subplot(1,2,cond); hold on
            a = plot(phiBinMid,nanmean(y),'color',cmap(5,:),'LineWidth',2);
            if cond == 1; a.LineStyle = '--'; end
    end
    
end




