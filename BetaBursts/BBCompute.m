function [BB flag] = BBCompute(R,surflag,sub,side,breg)

BB = []; flag = 0;
for cond = 1:length(R.condname)
    if cond<3
        if ~strncmp(R.subname{sub},'OXSH_D12',3)
            load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
        else
            load([R.datapathr R.subname{sub} '\ftdata\cleaned\V6_sources_clean_ROI_' R.condname{1} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'vc_clean')
        end
        else
        % %                         surflag = 1;
    end
    if  vc_clean.specanaly.flag == 1
        flag = flag +1 ;
        break
    end
    if cond == 1
        BB = getBurstSyncFrqs(R,vc_clean,BB);
    end
    BB = compute_BetaBursts(R,BB,vc_clean,cond,surflag);
    % BB.APrc{cond}
    BB.fsamp = vc_clean.fsample;
    
end