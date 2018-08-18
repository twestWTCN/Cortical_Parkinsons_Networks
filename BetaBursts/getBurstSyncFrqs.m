function BB = getBurstSyncFrqs(R,data,BB)
if nargin<3
    BB = [];
end
bind = find(data.specanaly.frq >= R.bandef(2,1) &  data.specanaly.frq <= R.bandef(3,2));
[pm pi] = max(data.specanaly.normpow(2,bind));
[cm ci] = max(data.specanaly.coh(1,bind));
BB.powfrq = data.specanaly.frq(bind(pi));
BB.cohfrq = data.specanaly.frq(bind(ci));