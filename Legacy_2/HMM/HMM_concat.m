function HMM_concat(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end
%%%
close all
for sub = [1:3 5:length(R.subname)]
    for breg = 2:length(R.bregname)
        for cond = 1:length(R.condname)
            clear jhem
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                x = vc_clean.trial{1};
                %                 y = (x-mean(x,2))./std(x,1,2);
                jhem{side} = x;
                
                data{sub,side,cond} = horzcat(jhem{:})';
                T{sub,side,cond} = size(horzcat(jhem{:})',1);
            end
        end
    end
end
data = reshape(data,[],1);
T = reshape(T,[],1);

% ! shutdown /h
mkdir([R.datapathr 'HMMdata\HMMdataconcat'])
save([R.datapathr 'HMMdata\HMMdataconcat'],'data','T')
a = 1;