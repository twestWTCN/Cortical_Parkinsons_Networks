function BB = getBurstSyncFrqs(R,data,BB)
if nargin<3
    BB = [];
end
bind = find(data.specanaly.frq >= R.bandef(2,1) &  data.specanaly.frq <= R.bandef(2,2)); % Find the beta range (B1 and B2)
[pm pi] = max(data.specanaly.normpow(2,bind)); % Find the peak in the STN power spectra
bind = find(data.specanaly.frq >= R.bandef(2,1) &  data.specanaly.frq <= R.bandef(3,2)); % Find the beta range (B1 and B2)
[cm ci] = max(data.specanaly.coh(1,bind));     % Find the peak in the STN/SMA coh spectra
BB.powfrq = data.specanaly.frq(bind(pi));
BB.cohfrq = data.specanaly.frq(bind(ci));