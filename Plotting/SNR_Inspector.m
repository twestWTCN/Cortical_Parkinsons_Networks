function SNR_Inspector(R,timev,SNR_sw,SNR_eps_z,breg)
for i = 1:numel(R.bandname)
    cmap = linspecer(3);
    ax(i) = subplot(1,numel(R.bandname),i);
    a(1) = plot(timev,SNR_sw(:,i),'color',cmap(i,:)); hold on;
    a(2) = plot(timev([1 end]),[SNR_eps_z(i) SNR_eps_z(i)],'color',cmap(i,:));
    xlabel('Time'); ylabel('dB'); ylim([-15 5]);% legend({'Signal',[R.PA.SNReps_prctile 'th Percentile']})
    if i == 1
    title([R.bregname{breg} ' ' R.bandname{R.bregband{breg}} ' band SNR'])
    elseif i == 2
    title([R.titular{2} ' ' R.bandname{R.bregband{breg}} ' band SNR'])
    elseif i == 3
    title([R.titular{2} ' ' R.bandname{strncmp(R.bandname,'Low Beta',6)} ' band SNR'])
    end
    grid on 
end
linkaxes(ax,'xy');
set(gcf,'Position',[351         771        1520         327]); shg