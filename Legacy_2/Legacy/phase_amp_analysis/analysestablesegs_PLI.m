function [phi_dist amp_dist seg_ddt segL_ddt consecSegs H] = analysestablesegs_PLI(qstable,tseries,refseries,period,mwid,fsamp,timevec)
amp = refseries;
consecSegs = SplitVec(qstable,'consecutive');
% lengths
segL_ddt = cellfun('length',consecSegs);
seglist = find(segL_ddt>(period));

seg_ddt = [consecSegs{segL_ddt>(period)}];
segL_ddt = segL_ddt(segL_ddt>(period))/fsamp;
% phase angles of segments
phi_dist = [];
amp_dist = [];
H = [];
for j = 1:numel(seglist)
        phi_dist(j) = circ_mean(wrapToPi(tseries(consecSegs{seglist(j)})));
%         phi_2_dist(j) = circ_mean(wrapTo2Pi(tseries(consecSegs{j},2)));
%         phi_1_dist(j) = circ_mean(wrapTo2Pi(tseries(consecSegs{j},1)));
        amp_dist(1,j) = (median(amp(consecSegs{seglist(j)},1))-median(amp(:,1)))/median(amp(:,1))*100;
        amp_dist(2,j) = (median(amp(consecSegs{seglist(j)},2))-median(amp(:,2)))/median(amp(:,2))*100;
        amp_dist(3,j) = (median(amp(consecSegs{seglist(j)},3))-median(amp(:,3)))/median(amp(:,3))*100;
        % Amp correlations
        x1 = amp(consecSegs{seglist(j)},1); x2 = amp(consecSegs{seglist(j)},2);
        [r,i] = my_xcorr(x1,x2,-256:256); %[r,i] = xcorr(x1,x2,fsamp,'unbiased');
         [r,ii] = max(r);
%         subplot(1,2,1); plot(1:numel(x1),x1,1:numel(x2),x2);
%         subplot(1,2,2); plot(x1+x2);
        H(1,j) = r;
        H(2,j) = i(ii)/fsamp;
%         H(3,j) = mean(amp(consecSegs{seglist(j)},3))/mean(diff(amp(consecSegs{seglist(j)},3)));
%         H(1,j) = mean(amp(consecSegs{seglist(j)},1))/mean(diff(amp(consecSegs{seglist(j)},1)));
%         H(2,j) = mean(amp(consecSegs{seglist(j)},2))/mean(diff(amp(consecSegs{seglist(j)},2)));
%         H(3,j) = mean(amp(consecSegs{seglist(j)},3))/mean(diff(amp(consecSegs{seglist(j)},3)));
% ####        amp_dist(3,j) = max(amp(consecSegs{j},3))/mean(amp(:,3));
end
% a= 1;
% phi_dist(isnan(phi_dist)) = [];
% amp_dist(isnan(amp_dist)) = [];
