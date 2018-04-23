clear
consecSegs{1} = 1:24;
consecSegs{2} = 25:144;
consecSegs{3} = 255:286;
consecSegs{4} = 299:301;
consecSegs{5} = 303:334;
consecSegs{6} = 340:382;
consecSegs{7} = 425:482;
% seglistCon{1} = [consecSegs{1} consecSegs{2}] 144
% seglistCon{2} = [consecSegs{3}] 32
% seglistCon{3} = [consecSegs{4} consecSegs{5} consecSegs{6}] 78
% seglistCon{4} = [consecSegs{7}]58
minsep = 8;
% j = 1;
% seglistCon{1} = consecSegs{1};
% i = 1;j = 2
% while j<numel(consecSegs)
%     if (consecSegs{j}(1)-seglistCon{i}(end))<minsep
%          seglistCon{i} = [seglistCon{i} consecSegs{j}]
%          j = j+2; i = i+1;
%     else
%         seglistCon{i+1} = consecSegs{j};
%         j = j+1; i = i+1;
%     end
% end

for i = 1:(numel(consecSegs)-1)
    gap(i) = consecSegs{i+1}(1) -consecSegs{i}(end)
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
   
        
        

% for j = 1:numel(consecSegs)-1
%     count = count+1;
%     if (consecSegs{j+1}(1)-consecSegs{j}(end))<minsep
%          seglistCon{count} = [consecSegs{j} consecSegs{j+1}]
%          j = j+1; % skip the next
%     else
%          seglistCon{count} = consecSegs{j}
%     end
%
%
% end
