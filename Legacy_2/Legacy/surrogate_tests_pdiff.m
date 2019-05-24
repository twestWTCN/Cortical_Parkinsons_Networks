clear; close all
load('Xdata.mat')
% load('surdata.mat')
Odata= Xdata;
% [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,37,16,0.5);
% qstable = find(abs(dphi_12_dt')<0.01); % The points that are below threshold
% consecSegs = SplitVec(qstable,'consecutive');
% % lengths
% segL_ddt = cellfun('length',consecSegs);
% histogram(log(segL_ddt/256),'FaceColor','b'); hold on

surdata(1,:) = Xdata.trial{1}(1,:);
surdata(2,:) = Xdata.trial{1}(2,:);


surfft1 = fft(surdata(1,:));
Amp1 = abs(surfft1);
Phi1 = angle(surfft1);

surfft2 = fft(surdata(2,:));
Amp2 = abs(surfft2);
Phi2 = angle(surfft2);


N = 6;
eps = logspace(-3,-1.5,N);
% eps(1) = 0.002;
% eps = 0.38;

WinN = [0.5 1 2.5]; %[0.05 0.1 0.25 0.5 2.5 3 5];
bwidlist = [0.5 1 2.5];
frqlist = 24:0.5:36;

for i = 1:N
    for bwid = 1:numel(bwidlist)
        for rep = 1:3
            parfor frqn = 1:numel(frqlist)
                %         figure
                [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,frqlist(frqn),16,bwidlist(bwid),1024);
                figure(1)
                histogram(abs(dphi_12_dt),logspace(-6,1.5,25),'FaceColor','b','Normalization','probability'); hold on; set(gca,'xscale','log'); grid on; ylabel('P(X)')
                qstable = find(abs(dphi_12_dt)<eps(i)); % The points that are below threshold
                consecSegs = SplitVec(qstable,'consecutive'); PLVO = abs(dphi_12_dt); real = mean(abs(dphi_12_dt));
                xlabel('Abs dphi_dt')
                % lengths
                segL_ddt = cellfun('length',consecSegs); segLVO = segL_ddt/1024;
                realL = mean(segL_ddt/1024);
                figure(2)
                histogram(segL_ddt/1024,logspace(-4.5,1.5,25),'FaceColor','b','Normalization','probability'); hold on; set(gca,'xscale','log'); grid on; ylabel('P(X)')
                xlabel('Segment Length (s)')
                figure(3)
                plot(Odata.time{1}(1024+1:end-1024),abs(dphi_12_dt),'b'); hold on; plot(Odata.time{1}(2:end),repmat(eps(i),1,length(Odata.time{1}(2:end))))
                plot(Odata.time{1}(1024+1:end-1024),repmat(mean(abs(dphi_12_dt)),1,length(dphi_12_dt)),'b--')
                
                figure(4)
                plot(betaS.time{1},betaS.trial{1}(1,:)); hold on
                
                Rdata = [];
                Rdata = Odata;
                Phishuff = Phi1(randperm(length(Phi1)));
                phasescm1 = abs(ifft(Amp1.*exp(sqrt(-1)*Phishuff))); %randn(1,length(Phishuff)); %
                Phishuff = Phi2(randperm(length(Phi2)));
                phasescm2 = abs(ifft(Amp2.*exp(sqrt(-1)*Phishuff)));
                
                Rdata.trial{1} = [phasescm1; phasescm2];
                
                [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frqlist(frqn),16,bwidlist(bwid),1024);
                figure(4)
                plot(Odata.time{1},betaS.trial{1}(1,:)); hold on

                figure(1)
                histogram(abs(dphi_12_dt),logspace(-6,1.5,25),'FaceColor','r','Normalization','probability'); hold on; set(gca,'xscale','log'); grid on; ylabel('P(X)')
                legend({'Real','Surr'})
                qstable = find(abs(dphi_12_dt)<eps(i)); % The points that are below threshold
                consecSegs = SplitVec(qstable,'consecutive'); PLV = abs(dphi_12_dt); sur = mean(abs(dphi_12_dt)); 
                % lengths
                segL_ddt = cellfun('length',consecSegs); segLV = segL_ddt/1024;
                surL = mean(segL_ddt/1024);
                figure(2)
                histogram(segL_ddt/1024,logspace(-4.5,1.5,25),'FaceColor','r','Normalization','probability'); hold on; set(gca,'xscale','log')
                legend({'Real','Surr'})
                figure(3)
                plot(Odata.time{1}(1024+1:end-1024),abs(dphi_12_dt),'r'); hold on; plot(Odata.time{1}(2:end),repmat(eps(i),1,length(Odata.time{1}(2:end))))
                plot(Odata.time{1}(1024+1:end-1024),repmat(mean(abs(dphi_12_dt)),1,length(dphi_12_dt)),'r--')
                
                epsi(frqn) = prctile(PLV,25);
                thresh(frqn) = prctile(segLV,75);
                [h pval_PLV(frqn)] = ttest2(log10(PLV),log10(PLVO));
                diffmeans(frqn) = sur-real;
                diffsegL(frqn) = surL-realL;
                [h pval_segL(frqn)] = ttest2(log10(segLV),log10(segLVO));
                close all
                disp([frqn])
            end
            epsi_save(:,rep,bwid,i) = epsi;
            thresh_save(:,rep,bwid,i) =thresh;
            pval_PLV_save(:,rep,bwid,i) =pval_PLV;
            pval_segL_save(:,rep,bwid,i) =pval_segL;
            diffmeans_save(:,rep,bwid,i) =diffmeans
            diffsegL_save(:,rep,bwid,i) =diffsegL
            disp([bwid i])
        end
    end
end
N = 6;
close all
for i= 1:N %numel(WinN)
    figure(1); %set(gcf,'Position',[1497          74         271         920])
    subplot(N,1,i)
    plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(pval_segL_save(:,:,:,i),2)),'LineWidth',4);% ylim([-0.01 0])
    legend({sprintf('%.1f Hz',bwidlist(1)),sprintf('%.1f Hz',bwidlist(2)),sprintf('%.1f Hz',bwidlist(3))})
    ylabel('SRP Length ttest P-Value'); xlabel('Freq'); title(['Eps ' num2str(eps(i))]); grid on;% set(gca,'yscale','log')
    figure(2);% set(gcf,'Position',[1497          74         271         920])
    subplot(N,1,i)
    plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(pval_PLV_save(:,:,:,i),2)),'LineWidth',4 );% ylim([-0.1 0])
    legend({sprintf('%.1f Hz',bwidlist(1)),sprintf('%.1f Hz',bwidlist(2)),sprintf('%.1f Hz',bwidlist(3))})
    ylabel('|dRP| ttest P-Value'); xlabel('Freq'); title(['Eps ' num2str(eps(i))]); grid on;% set(gca,'yscale','log')
    
    figure(3); %set(gcf,'Position',[1497          74         271         920])
    subplot(N,1,i)
    plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(epsi_save(:,:,:,i),2)),'LineWidth',4 ); %ylim([0.1 1])
    legend({sprintf('%.1f Hz',bwidlist(1)),sprintf('%.1f Hz',bwidlist(2)),sprintf('%.1f Hz',bwidlist(3))})
    ylabel('25% Prctile |dRP| rad.s^-1'); xlabel('Freq'); title(['Eps ' num2str(eps(i))]); grid on
    
    figure(4); %set(gcf,'Position',[1497          74         271         920])
    subplot(N,1,i)
    plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(thresh_save(:,:,:,i),2)),'LineWidth',4 ); %ylim([10 90])
    legend({sprintf('%.1f Hz',bwidlist(1)),sprintf('%.1f Hz',bwidlist(2)),sprintf('%.1f Hz',bwidlist(3))})
    ylabel('75% Prctile SRP Length (s)'); xlabel('Freq'); title(['Eps ' num2str(eps(i))]); grid on
    
end
savefigure_v2([R.datapathr '\results\seganalysis\surrogates\'],['PLI_surrogate'],[],[],[]); % close all

R.PA.bwid = 1;
R.PA.eps = 0.003;
R.PA.slidingwindow = 2.5;
R.PA.PLVeps =  0.35;
R.PA.mwid = 30;
R.PA.WinOver = 0.95;
R.PA.stn_lb_frq = 14;
