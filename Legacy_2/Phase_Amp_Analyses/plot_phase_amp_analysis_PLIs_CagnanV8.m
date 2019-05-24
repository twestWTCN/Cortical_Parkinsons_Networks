function plot_phase_amp_analysis_PLIs_CagnanV8(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
close all
analynames = {'Segment Length','CTX High Beta Amp','STN High Beta Amp','STN Low Beta Amp','STN/CTX High Beta Amp Correlation','Causal Density'};
% WIP
% Percentage Change over what? Baseline or baseline of segments?
QX = 8; % Bin Size
% sub 3 side 1
for breg = 2:length(R.bregname)
    cmapint = linspecer(5);
    cmap = linspecer(5);
    cmapint(5,:) = cmap(5,:);
    ampBinGroup = []; segBinGroup = []; phipeakGroup = []; RcoeffGroup = [];
    for sub = 1:length(R.subname)
        cmapint = cmapint*0.95;
        for side = 1:2
            for cond = 1:length(R.condname)
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                if vc_clean.specanaly.flag ~= 1 % Check significant coherences                        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\' idd '_ROIvoxel_phaseamp_PLI_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '.mat'])
                    phi = vc_clean.PA.pA_pli_dist_save';
                    relativePhiCol = phi'; %wrapToPi(phi-circ_mean(phi(~isnan(phi'))) )'; %circ_mean(pA_pli_dist_save{2,nrOFF}(amp_pli_dist_save{2,nrOFF}(3,:)>95)')); %wrapToPi(pA_pli_dist_save{cond,nr}); %-circ_mean([pA_pli_dist_save{2,nrOFF}]'));
                    segLCol =  vc_clean.PA.segL_pli_dist_save; %((segL_pli_dist_save{cond,nr} - mean([segL_pli_dist_save{cond,nr}],2))./mean([segL_pli_dist_save{cond,nr}],2)  )*100;
                    ampSegCol = vc_clean.PA.amp_pli_dist_save; %((amp_pli_dist_save{cond,nr}  - mean([amp_pli_dist_save{cond,nr}],2))./mean([amp_pli_dist_save{cond,nr}],2) )*100;
                    %                     HdistSegCol = vc_clean.PA.H_dist_save(1,:);
                    tendtot = vc_clean.PA.timevec{1}(end)-vc_clean.PA.timevec{1}(1);
                    
                    ampdata{cond,side,sub} = ampSegCol;
                    segdata{cond,side,sub} = segLCol;
                    rpdata{cond,side,sub} = relativePhiCol;


                    ylimlistS{1}= {{[-25 350];[-25 250];[-25 250];[-25 250]},{[-25 100];[-25 300];[-25 100];[-25 300]}};
                    ylimlistS{2}={{[0 1]},{[0 1.5]},{[0 1]}};
                    ylimlistS{3}= {{[-25 350];[-25 250];[-25 250];[-25 250]},{[-100 300];[-100 400];[-100 400];[-100 400]}};
                    panlist = [1 3 5 7 9; 2 4 6 8 10];
                    phiBin = linspace(-pi,pi,QX);
                    obs = {[R.bregname{breg} ' ' R.bandinits{R.bregband{breg}}],['STN ' R.bandinits{R.bregband{breg}}],...
                         [R.bregname{breg} ' '  R.bandinits{2}],['STN ' R.bandinits{2}],[R.bregname{breg} ' Seg. Length']};
                    for i = 1:4
                        titleR{1,i} = ['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{i} ' Power'];
                        titleR{2,i} =  ['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{i} ' Power'];
                        titleR{3,i} =['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' R.bandinits{R.bregband{breg}} ' Frame Length'];
                    end
                    titleR{4,4} = ['STN-' R.bregname{breg} ' ' R.bandinits{R.bregband{breg}} ' Phase vs ' obs{5}];
                    [ampBinGroup segBinGroup phipeakGroup RcoeffGroup] =  plot_phase_amp_analysis_generic(R,ylimlistS,obs,breg,cond,side,sub,phiBin,...
                        cmapint,segLCol,relativePhiCol,ampSegCol,ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup,...
                        panlist,titleR);
                ampsurdata{cond,side,sub} = makeAmpRPhiSurr(R,ampSegCol,relativePhiCol,phiBin);

                end
            end
                    close all
        end
    end
    plot_phase_amp_analysis_generic_group(R,panlist,obs,phiBin,ampBinGroup,segBinGroup,phipeakGroup,RcoeffGroup,ampdata,segdata,rpdata,ampsurdata)
    close all
end

