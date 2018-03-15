clear; %close all
load('Xdata.mat')
% load('surdata.mat')
Odata= Xdata;
% [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,37,16,0.5);
% qstable = find(abs(dphi_12_dt')<0.01); % The points that are below threshold
% consecSegs = SplitVec(qstable,'consecutive');
% % lengths
% segL_ddt = cellfun('length',consecSegs);
% histogram(log(segL_ddt/256),'FaceColor','b'); hold on

surdata(1,:) = Xdata.trial{1}(:,:);
surdata(2,:) = Xdata.trial{1}(2,:);


surfft1 = fft(surdata(1,:));
Amp1 = abs(surfft1);
Phi1 = angle(surfft1);

surfft2 = fft(surdata(2,:));
Amp2 = abs(surfft2);
Phi2 = angle(surfft2);


N = 6;
eps = logspace(-3,-2,N);
eps = 0.38;

WinN = [0.5 1 2.5]; %[0.05 0.1 0.25 0.5 2.5 3 5];
bwidlist = [0.5 1 2.5];
frqlist = 24:0.5:36;

for i = 1:numel(WinN)
    for bwid = 1:numel(bwidlist)
        for rep = 1:3
        for frqn = 1:numel(frqlist)
            %         figure
            [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Odata,frqlist(frqn),16,bwidlist(bwid));
            WinSize =WinN(i)*1024;
            [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,0.95);
            real = mean(PLV);
            SW_sampr = max(diff(PLV_tvec))/1024;
            %         figure(1)
            %         histogram(PLV,'FaceColor','b','Normalization','pdf'); hold on
            qstable = find(PLV>eps); % The points that are below threshold
            consecSegs = SplitVec(qstable,'consecutive'); PLVO = PLV;
            % lengths
            segL_ddt = cellfun('length',consecSegs);
            realL = mean(segL_ddt/1024);
            %         figure(2)
            %         histogram(segL_ddt/SW_sampr,'FaceColor','b','Normalization','pdf'); hold on
            %         figure(3)
            %         plot(PLV_tvec/1024,PLV,'b'); hold on; plot(PLV_tvec/1024,repmat(eps,1,numel(PLV_tvec)))
            %         plot(PLV_tvec/1024,repmat(mean(PLV),1,length(PLV)),'b--')
            
            Phishuff = Phi1(randperm(length(Phi1)));
            phasescm1 = randn(1,length(Phishuff)); %real(ifft(Amp1.*exp(sqrt(-1)*Phishuff)));
            Phishuff = Phi2(randperm(length(Phi2)));
            phasescm2 = randn(1,length(Phishuff)); %real(ifft(Amp2.*exp(sqrt(-1)*Phishuff)));
            Rdata = [];
            Rdata.trial{1} = [phasescm1; phasescm2];
            
            [amp phi dphi_12 dphi_12_dt betaS] = comp_instant_angle_phase(Rdata,frqlist(frqn),16,bwid);
            WinSize = WinN(i)*1024;
            [PLV PLV_tvec] = slidingwindowPLV(WinSize,phi,0.95);
            %         figure(1)
            %         histogram(PLV,'FaceColor','r','Normalization','pdf'); hold on
            
            qstable = find(PLV>eps); % The points that are below threshold
            consecSegs = SplitVec(qstable,'consecutive');
            segL_ddt = cellfun('length',consecSegs);
            surL = mean(segL_ddt/1024);
            %         figure(2)
            %         histogram(segL_ddt/1024,'FaceColor','r','Normalization','pdf'); hold on
            %         figure(3)
            %         plot(PLV_tvec/1024,PLV,'r'); hold on; plot(PLV_tvec/1024,repmat(eps,1,numel(PLV_tvec)));
            %         plot(PLV_tvec/1024,repmat(mean(PLV),1,length(PLV)),'r--')
            sur = mean(PLV);
            
            epsi(frqn) = prctile(PLV,95);
            thresh(frqn) = prctile(segL_ddt,95);
            [h pval(frqn)] = ttest(PLV,PLVO);
            diffmeans(frqn) = sur-real;
            diffsegL(frqn) = surL-realL;
            close all
            disp([frqn])
        end
        epsi_save(:,rep,bwid,i) = epsi;
        thresh_save(:,rep,bwid,i) =thresh;
        pval_save(:,rep,bwid,i) =pval;
        diffmeans_save(:,rep,bwid,i) =diffmeans
        diffsegL_save(:,rep,bwid,i) =diffsegL
        disp([bwid i])
        end
    end
end

close all
for i= 1:numel(WinN)
 figure(1); set(gcf,'Position',[680.0000   58.5000  674.5000  919.5000]) 
subplot(3,1,i)
plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(diffsegL_save(:,:,:,i),2)) ); ylim([-0.01 0])
ylabel('Difference in SegL'); xlabel('Freq'); title(['WinSize ' num2str(WinN(i))])
 figure(2); set(gcf,'Position',[680.0000   58.5000  674.5000  919.5000]) 
subplot(3,1,i)
plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(diffmeans_save(:,:,:,i),2)) ); ylim([-0.1 0])
ylabel('Difference in Mean PLV'); xlabel('Freq'); title(['WinSize ' num2str(WinN(i))])

 figure(3); set(gcf,'Position',[680.0000   58.5000  674.5000  919.5000]) 
subplot(3,1,i)
plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(epsi_save(:,:,:,i),2)) ); ylim([0.1 1])
ylabel('95% Prctile PLV'); xlabel('Freq'); title(['WinSize ' num2str(WinN(i))])

 figure(4); set(gcf,'Position',[680.0000   58.5000  674.5000  919.5000]) 
subplot(3,1,i)
plot(repmat(frqlist,numel(WinN),1)',squeeze(mean(thresh_save(:,:,:,i),2)) ); ylim([10 90])
ylabel('95% Prctile SegL'); xlabel('Freq'); title(['WinSize ' num2str(WinN(i))])

end
savefigure_v2([R.datapathr '\results\seganalysis\surrogates\'],['PLI_surrogate'],[],[],[]); % close all

R.PA.bwid = 1;
R.PA.slidingwindow = 2.5;
R.PA.PLVeps =  0.35;
R.PA.mwid = 30;
R.PA.WinOver = 0.95;
R.PA.stn_lb_frq = 14;
