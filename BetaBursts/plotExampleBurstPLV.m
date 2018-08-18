function plotExampleBurstPLV(R,BB)
set(gcf,'defaultAxesColorOrder',[[0 0 0]; [0.2 0.2 0.6]]);
for i = 1:3
    if i==1
        subplot(2,1,1)
    else
        subplot(2,2,i+1)
    end
    yyaxis left
    a(1) = plot(BB.T,[BB.A{1}(2,:) BB.A{2}(2,:)]); hold on
    a(1).Color = [0 0 0]; ylim([0 125])
    a(3) = plot(BB.T([1 end]),[BB.epsAmp(2) BB.epsAmp(2)],'--');
    a(3).Color = [0 0 0];
    ylabel([R.bandinits{2} ' Amplitude'])
    % Plot Sync
    yyaxis right;
    a(2) = plot(BB.TSw,[BB.PLV{1}(3,:) BB.PLV{2}(3,:)]); hold on
    a(2).Color = [0.2 0.2 0.6]; ylim([0 1]);
    a(4) = plot(BB.T([1 end]),[BB.epsPLV(2) BB.epsPLV(2)],'--');
    a(4).Color = [0.2 0.2 0.6];
    yl = ylabel([R.bandinits{3} ' PLV']);
%     yl.Rotation = -90;
%     set(yl, 'Units', 'Normalized', 'Position', [1.12, 0.5, 0]);
    xlabel('Time (s)');
    
    % Plot the OFF Condition Marker
    yyaxis left
    plot(BB.DTvec{1},repmat(120,size(BB.DTvec{1})),'k-','LineWidth',3,'Color',R.condcmap(1,:))
    box off
    if i==1;     title('Example Wavelet Analysis');
        legend(a(1:2),{[R.bandinits{2} ' Amplitude'],[R.bandinits{3} ' PLV']},'Location','NorthEast'); end
    if i==2; for p=1:4; a(p).LineWidth = 1.5; end; xlim([50 60]); title(R.condname{i-1}); end
    if i==3; for p=1:4; a(p).LineWidth = 1.5; end; xlim([250 260]); title(R.condname{i-1}); end
end
