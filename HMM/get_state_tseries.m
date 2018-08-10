function Ts = get_state_tseries(data,S,trange)
clear stTs
X = vertcat(data{:});
for st = 1:size(S,2)
    i = 0;
    stTs = {};
    for tr = trange
        for visit = 1:size(S{tr,st},1)
            i = i+1;
            stTs{i} = X(S{tr,st}(visit,1):S{tr,st}(visit,2),:);
        end
    end
    Y = vertcat(stTs{:});
    Ts{st} = Y;
end