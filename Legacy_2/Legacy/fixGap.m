function [seglistCon] = fixGap(consecSegs,minsep)

for i = 1:(numel(consecSegs)-1)
    gap(i) = consecSegs{i+1}(1) -consecSegs{i}(end);
end
mingap = gap<minsep;
i = 1;
seglistCon{1} = consecSegs{1};
for j = 1:length(mingap)
    if mingap(j) == 1
        seglistCon{i} = [seglistCon{i} consecSegs{j+1}];
    else
         i = i+1; j = j+1;
        seglistCon{i} = consecSegs{j};
    end
end