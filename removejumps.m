function [X_clean] = removejumps(X,thresh,plotop)
X_clean = X;
spt = 10; % width of cleaning region
bksmp1 = 5; % back sampling
XThresh = find(abs(X)>thresh*std(X));
consecSegs = SplitVec(XThresh,'consecutive');
T = 1:numel(X);
if plotop; plot(T,X); hold on; end
for seg = consecSegs
    if seg{1}(end)+(spt+bksmp1)>length(X)
        bksmp = length(X)-seg{1}(end)-spt;
    else
        bksmp = bksmp1;
    end
    if length(X)-seg{1}(end)<bksmp
        X([seg{1}(1)-spt:seg{1}(end)]) = NaN;
    elseif seg{1}(1)<=bksmp
        X([seg{1}(1):seg{1}(end)+spt]) = NaN;
    else
        try
            X([seg{1}(1)-spt:seg{1}(end)+spt]) = NaN;
        catch
            X([seg{1}(1)-(seg{1}(1)-1):seg{1}(end)+spt]) = NaN;
        end
    end
end
consecY = SplitVec(find(isnan(X)),'consecutive');

for seg = consecY
    if seg{1}(1)-bksmp<1
        while seg{1}(1)-bksmp<1;        bksmp =  bksmp-1;            end
    end
    if seg{1}(end)+bksmp>numel(X)
        while seg{1}(end)+bksmp>numel(X);   bksmp = bksmp-1;    end
    end
    if seg{1}(end)+(spt+bksmp1)>length(X)
        bksmp = length(X)-seg{1}(end)-spt;
    else
        bksmp = bksmp1;
    end
    if length(X)-seg{1}(end)<bksmp
        Xrep = interp1(T([seg{1}(1)-(spt+bksmp):seg{1}(1)-1]),...
            X([seg{1}(1)-(spt+bksmp):seg{1}(1)-1]),...
            T([seg{1}(1)-1:seg{1}(end)]),'spline');
        Tind = T([seg{1}(1)-1:seg{1}(end)]);
    elseif seg{1}(end)>=length(X)
        Xrep = nanmean(X);
    elseif seg{1}(1)<bksmp
        Xrep = interp1(T([seg{1}(end)+1:seg{1}(end)+(spt+bksmp)]),...
            X([seg{1}(end)+1:seg{1}(end)+(spt+bksmp)]),...
            T([seg{1}(1):seg{1}(end)+1]),'spline');
         Tind = T([seg{1}(1):seg{1}(end)+1]);

    elseif length(T([seg{1}(1)-(bksmp-1):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+1]))<2
        Xrep = nanmean(X);
    else
        Xrep = interp1(T([seg{1}(1)-(bksmp-1):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+1]),...
            X([seg{1}(1)-(bksmp-1):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+1]),...
            T([seg{1}(1)-1:seg{1}(end)+1]),'spline');
        Tind = T([seg{1}(1)-1:seg{1}(end)+1]);
    end
    
    Xrep = Xrep - nanmean(X);
    flag = 0; Xlen = numel(Xrep);
    while flag ==0
        rndst = ceil(rand.*T(end));
        if rndst<1; rndst = 1; end
        try
            C = X(rndst:rndst+Xlen-1);
        catch
            C = X(rndst:end);
            disp('Exceeds limits, filling to end of signal!!')
        end
        if ~any(isnan(C))
            Xrep = sign(nanmean(Xrep)).*C;
            flag = 1;
        else
            flag = 0;
        end
    end
    Tind = Tind(1:length(Xrep));
    X_clean(Tind) = Xrep;
    if plotop; plot(T(Tind),Xrep,'r'); end
end
if plotop; plot(T,X_clean,'g'); end
