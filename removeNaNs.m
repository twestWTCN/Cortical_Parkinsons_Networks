function [X_clean flag] = removeNaNs(X,plotop)
flag = 0;
bksmp = 24; % back sampling
ARsmp = 512;
XNaN = find(isnan(X));
if numel(XNaN)>(numel(X)*0.4)
    flag = 1;
    X_clean = X;
    return
end
consecSegs = SplitVec(XNaN,'consecutive');
T = 1:numel(X);
if plotop; plot(T,X); hold on; end
for seg = consecSegs
    %% Interpolation
    try
        Xrep = interp1(T([seg{1}(1)-(bksmp):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+(bksmp)]),...
            X([seg{1}(1)-(bksmp):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+(bksmp)]),...
            T([seg{1}(1)-1:seg{1}(end)+1]),'linear');
    catch
        if seg{1}(1)<2
            disp('No interpolation possible, assuming linearity')
            Xrep = linspace(mean(seg{1}),seg{1}(end),numel(seg{1}));
        elseif seg{1}(end)>length(X)
            Xrep = linspace(seg{1}(1),X(end),numel(seg{1}(1):X(end)));
        else
            if seg{1}(1)-bksmp<1
                while seg{1}(1)-bksmp<1;        bksmp =  bksmp-1;            end
            end
            if seg{1}(end)+bksmp>numel(X)
                while seg{1}(end)+bksmp>numel(X);   bksmp = bksmp-1;    end
            end
            Xrep = interp1(T([seg{1}(1)-(bksmp):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+(bksmp)]),...
                X([seg{1}(1)-(bksmp):seg{1}(1)-1, seg{1}(end)+1:seg{1}(end)+(bksmp)]),...
                T([seg{1}(1)-1:seg{1}(end)+1]),'linear');
        end
    end
    %% AR Model
    if seg{1}(1)<2 || seg{1}(end)>length(X)
        a = rand(size(Xrep));
        X = Xrep.*a;
    else
        m = ar(X([seg{1}(1)-(ARsmp):seg{1}(1)-1]),32);
        a = sim(m,rand(1,ceil(numel(Xrep)*1.1))');
        a = a(end-length(seg{1})-1:end)'; a = a-mean(a);
        Xrep = Xrep+(sign(mean(Xrep))*a); %sign(mean(Xrep))*sqrt(abs(Xrep)); %.*X(rndst:rndst+numel(Xrep)-1);
        X([seg{1}(1)-1:seg{1}(end)+1]) = Xrep;
    end
    if plotop;  plot(T([seg{1}(1)-1:seg{1}(end)+1]),Xrep,'r'); end
    
end
X_clean = X;
if plotop; plot(T,X_clean,'g'); end
close all