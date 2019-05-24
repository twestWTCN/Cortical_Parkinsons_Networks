function plotAmpBinDurBar(R,BB)
pvec = pvec_bin_TTest(Seg_binAmpDataLOC);
barplot160818(R,BB.binAmp,squeeze(BB.Seg_binAmp(:,:,:)),pvec,0,1)
title('Burst Duration by Amplitude')
xlabel('Mean Amplitude'); ylabel('Wghtd. Duration (ms)'); ylim([0 200]); % ylim([BB.range.Amp(1) BB.range.Amp(end)]);
