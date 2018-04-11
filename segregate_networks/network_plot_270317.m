clear; close all
R = makeHeader_SubCort_Cort_Networks_SFLAP();
% The normalization to zscore cant be used as a threshold as the smallest
% value is added to rescale. Consider saving that value. 
% Find the correct sampling frequencies and ensure allignment of time
% seires.
for sub = 1:6; %numel(R.subname)
    for cond = 1:2
        for side = 1:2
            p = 0; clear y2
            for band = [1 3]; %1:numel(R.bandname)
                p = p+1;
                load([R.datapathr R.subname{sub} '\ftdata\networks\threshnetwork_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '_' R.condname{cond}],'threshnetwork')
                %                 a(1) =subplot(3,1,1)
                %                 plot(threshnetwork.Xdata.time{1},log10(threshnetwork.amp(:,1))); hold on
                %                 a(2)= subplot(3,1,2)
                %                 plot(threshnetwork.Xdata.time{1},log10(threshnetwork.amp(:,2))); hold on
                
                a(1) =subplot(3,1,2)
                b1(p) = plot(threshnetwork.PLV_tvec,threshnetwork.PLV); hold on
                y1{band} = threshnetwork.PLV;
                a(2) =subplot(3,1,1)
                tvec = linspace(0,size(threshnetwork.amp,1)/1024,size(threshnetwork.amp,1));
                y2{band} =normaliseV(threshnetwork.amp(:,1)');
                y2{band} = y2{band} + abs(min(y2{band}));
                y2{band} = movmean(y2{band},2*1024);
                b2(p) = plot(tvec,y2{band}); hold on
                %                 linkaxes(a,'x'); %xlim([80 82])
                
                %                 threshnetwork.PLV = PLV;
                %                 threshnetwork.PLV_tvec = PLV_tvec;
                %                 threshnetwork.stn_lb_frq = stn_lb_frq;
                %                 threshnetwork.maxfrq = maxfrq;
                %                 threshnetwork.betaS = betaS;
                %                 threshnetwork.amp = amp;
                %                 threshnetwork.Xdata = Xdata;
                %                 threshnetwork.consecSegs = consecSegs;
                %                 mkdir([R.datapathr R.subname{sub} '\ftdata\networks\'])
                %                 save([R.datapathr R.subname{sub} '\ftdata\networks\threshnetwork_' R.ipsicon '_' R.siden{side} '_' R.bandname{band} '_' R.condname{cond}],'threshnetwork')
            end
            clear tvecAB R_alpha_beta
            % Amp Env Corrcoef
            winLength =2;
            x1 = y2{1}; x2 = y2{3};
            sfrq = 1/median(diff(tvec));
            [slidex1,~] = slideWindow(x1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
            [slidex2,sind] = slideWindow(x2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
            z1 = (x1-mean(mean(slidex1,1)))./mean(std(slidex1,[],1));
            z2 = (x2-mean(mean(slidex2,1)))./mean(std(slidex2,[],1));
            for i = 1:size(slidex1,2)-1
                if median(z1(sind(:,i)))>10  || median(z2(sind(:,i)))>10
                    R_alpha_beta(i) = NaN;
                    tvecAB(i) = tvec(floor(median(sind(:,i))));
                else
                    sR = corrcoef(slidex1(:,i),slidex2(:,i));
                    R_alpha_beta(i) = sR(2);
                    tvecAB(i) = tvec(floor(median(sind(:,i))));
                end
            end
            R_alpha_bet_emp = R_alpha_beta;
            a(3) =subplot(3,1,3)
            b3(1) = plot(tvecAB,R_alpha_bet_emp,'b'); hold on
            
            clear R_alpha_beta tvecAB Rstore tstore
            parfor rep = 1:100
                R_alpha_beta = [];tvecAB = [];
                shuffX1 =   phaseperm(x1); %x1(randperm(length(x1)));
                shuffX2 =   phaseperm(x2); %x2(randperm(length(x2)));
                [slidex1,~] = slideWindow(shuffX1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
                [slidex2,sind] = slideWindow(shuffX2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
                
                for i = 1:size(slidex1,2)-1
                    sR = corrcoef(slidex1(:,i),slidex2(:,i));
                    R_alpha_beta(i) = sR(2);
                    tvecAB(i) = tvec(min(sind(:,i)));
                end
                Rstore(:,rep) = R_alpha_beta;
                tstore(:,rep) = tvecAB;
                rep
            end
            upperBndEnv = prctile(Rstore(:),95); upperBnd2 = prctile(Rstore(:),99);
            lowerBndEnv = prctile(Rstore(:),5); lowerBnd2 = prctile(Rstore(:),1);
            
% % %             upperBndEnv = 0.8766; tvecAB = tvecAB';
% % %             lowerBndEnv = -0.7857;
            
            b3(2) = plot(tstore(:,1),repmat(upperBndEnv,1,size(tstore,1)),'b--'); hold on
            b3(3) = plot(tstore(:,1),repmat(lowerBndEnv,1,size(tstore,1)),'b--');
            
            AmpEnvCorr(:,1) = sum(R_alpha_bet_emp>upperBndEnv).*median(diff(tstore(:,1)));
            AmpEnvCorr(:,2) = sum(R_alpha_bet_emp<lowerBndEnv).*median(diff(tstore(:,1)));
            AmpEnvCorr_save{sub,side,cond} = AmpEnvCorr;
            
            clear R_alpha_beta tvecAB
            %% PLV Corrcoef
            clear tvecAB R_alpha_beta
            winLength =2;
            x1 = y1{1}; x2 = y1{3};
            sfrq = 1/median(diff(threshnetwork.PLV_tvec));
            [slidex1,~] = slideWindow(x1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
            [slidex2,sind] = slideWindow(x2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
            for i = 1:size(slidex1,2)-1
                if median(z1(sind(:,i)))>10  || median(z2(sind(:,i)))>10
                    R_alpha_beta(i) = NaN;
                    tvecAB(i) = tvec(median(sind(:,i)));
                else
                    sR = corrcoef(slidex1(:,i),slidex2(:,i));
                    R_alpha_beta(i) = sR(2);
                    tvecAB(i) = threshnetwork.PLV_tvec(median(sind(:,i)));
                end
            end
            tvecAB = tvecAB + threshnetwork.PLV_tvec(1)
            a(3) =subplot(3,1,3)
            b3(4) = plot(tvecAB,R_alpha_beta,'k'); hold on
            RPLV_alpha_bet_emp = R_alpha_beta;
            
            
            clear R_alpha_beta tvecAB Rstore tstore
            parfor rep = 1:100
                R_alpha_beta = []; tvecAB = [];
                shuffX1 =   phaseperm(x1); %x1(randperm(length(x1)));
                shuffX2 =   phaseperm(x2); %x2(randperm(length(x2)));
                [slidex1,~] = slideWindow(shuffX1, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
                [slidex2,sind] = slideWindow(shuffX2, floor(winLength*sfrq), floor(winLength*sfrq*0.95));
            
                for i = 1:size(slidex1,2)-1
                    sR = corrcoef(slidex1(:,i),slidex2(:,i));
                    R_alpha_beta(i) = sR(2);
                    tvecAB(i) = threshnetwork.PLV_tvec(min(sind(:,i)));
                end
                Rstore(:,rep) = R_alpha_beta;
                tstore(:,rep) = tvecAB;                
                rep
            end
            upperBndPLV = prctile(Rstore(:),95); upperBnd2 = prctile(Rstore(:),99);
            lowerBndPLV = prctile(Rstore(:),5); lowerBnd2 = prctile(Rstore(:),1);
% % %             upperBndPLV = 0.7486; tvecAB = tvecAB';
% % %             lowerBndPLV = -0.7193;
            b3(5) = plot(tstore(:,1),repmat(upperBndPLV,1,size(tstore,1)),'k--'); hold on
            b3(6) = plot(tstore(:,1),repmat(lowerBndPLV,1,size(tstore,1)),'k--');
            % b2(4) = plot(tvecAB(:,1),repmat(upperBnd2,1,size(tvecAB,1)),'k--'); hold on
            % b2(5) = plot(tvecAB(:,1),repmat(lowerBnd2,1,size(tvecAB,1)),'k--');
            PLVCorr(:,1) = sum(RPLV_alpha_bet_emp>upperBndEnv).*median(diff(tstore(:,1)));
            PLVCorr(:,2) = sum(RPLV_alpha_bet_emp<lowerBndEnv).*median(diff(tstore(:,1)));
            PLVCorr_save{sub,side,cond} = PLVCorr;
            
            linkaxes(a,'x'); %xlim([80 82])
            
            subplot(3,1,2)
            xlabel('Frame Time (s)'); ylabel('Within-Frame PLV');
            legend(b1,{'Alpha STG Source','Beta Premotor Source'})
            subplot(3,1,1)
            xlabel('Frame Time (s)'); ylabel('Amplitude Envelope');
            legend(b2,{'Alpha STN/STG Network','PLV STN/Premotor Network'})
            subplot(3,1,3)
            xlabel('Frame Time (s)'); ylabel('Correlation Coeff.');
            legend(b3,{'Amp Env Correlation R','95th Percentile','5th Percentile','PLV Correlation R','95th Percentile','5th Percentile'}); %,'99th Percentile','1st Percentile'})
            
            set(gcf,'Position',[4141 -133 1302 821])
            close all
        end
    end
end

figure(10)
x1  = vertcat(AmpEnvCorr_save{:,:,1});
x2  = vertcat(AmpEnvCorr_save{:,:,2});
subplot(1,2,1)
barplot_N2(x1(:,1),x2(:,1),'','')
subplot(1,2,2)
barplot_N2(x1(:,2),x2(:,2),'','')

figure(20)
x1  = vertcat(PLVCorr_save{:,:,1});
x2  = vertcat(PLVCorr_save{:,:,2});
subplot(1,2,1)
barplot_N2(x1(:,1),x2(:,1),'','')
subplot(1,2,2)
barplot_N2(x1(:,2),x2(:,2),'','')


