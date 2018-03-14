R.PA.bwid = 0.5;
R.PA.mwid = 18;
R.PA.PLVeps = 0.3;

bwidLS = [0.5 1 1.5 2.5];
mwidLS = [8 12 16 20];
PLVepsLS = [0.1 0.2 0.3];
for i = 1:numel(PLVepsLS)
    for j = 1:numel(mwidLS)
        parfor k = 1:numel(bwidLS)
            Rtemp = R;
            Rtemp.PA.bwid = bwidLS(k);
            Rtemp.PA.mwid = mwidLS(j);
            Rtemp.PA.PLVeps = PLVepsLS(i);
            idd = sprintf('%.0f_%.f_%.f',i,j,k);
            compute_phase_amp_analysis_120318(Rtemp,idd)
            dwell = getHists_phase_amp_analysis_PLIs(Rtemp,idd)
            dwellsave{i,j,k} = dwell;
            disp(idd)
        end
        save([R.datapathr 'dwellstats'],'dwellsave')
    end
end
            
for i = 1:numel(PLVepsLS)
    for j = 1:numel(mwidLS)
        for k = 1:numel(bwidLS)
            p(i,j,k) = dwellsave{i,j,k}{3};
        end
    end
end