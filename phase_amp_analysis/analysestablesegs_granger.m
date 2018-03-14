function [gc_dist] = analysestablesegs_granger(qstable,tseries,refseries,period,mwid,fsamp)
amp = refseries;
consecSegs = SplitVec(qstable,'consecutive');
% lengths
segL_ddt = cellfun('length',consecSegs);
seglist = find(segL_ddt>(period));

seg_ddt = [consecSegs{segL_ddt>(period)}];
segL_ddt = segL_ddt(segL_ddt>(period))/fsamp;
% Granger causality of segments
progressStepSize = 1;
% parpool
ppm = ParforProgMon('Power Estimation: ', numel(seglist), progressStepSize, 800, 300);
gc_dist = zeros(9,numel(seglist));
parfor j = 1:numel(seglist)
        try
            spgc = compute_mvgc_pair_granger(tseries(consecSegs{seglist(j)},:),fsamp,[10 40],0);
            if isreal(spgc)
                gc_dist(:,j) = spgc(:);
            else
                error('Granger is Complex!!')
            end
        catch
            disp('Granger Failed!!')
            gc_dist(:,j) = NaN(9,1);
        end
        % ####        amp_dist(3,j) = max(amp(consecSegs{j},3))/mean(amp(:,3));
    ppm.increment();
end
ppm.delete()

