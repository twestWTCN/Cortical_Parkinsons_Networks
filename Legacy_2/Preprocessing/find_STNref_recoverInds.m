function find_STNref_recoverInds(R)
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
close all
for sub = 1:length(R.subname)
    for cond = 1:length(R.condname)
        for breg = 1:length(R.bregname)
            for side = 1:2
                load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
                chsellab{sub,cond,side} = vc_clean.label{2};
                
                
                disp([side sub cond])
                close all
            end
        end
    end
end

a