function CompileData(R)
if nargin<2
    R = makeHeader_SubCort_Cort_Networks();
end

%%%
% This script computes the instantaneous phase sync analyses. We partition
% the recording into stable segments using threshold, can compute Granger
% etc within frames. DFA is also computed PS/AE versions. Barplots for DFA
% at end of script.
%%%
for band = 3 %[1 3] 
    subi = 0;
    for sub = 1:numel(R.subname)
        load([R.datapathr 'results\spectral\screen.mat'])
        load([R.datapathr 'results\spectral\refidSave_' R.bandname{band} '.mat'])
        load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_' R.bandname{band}])
        for side = 1:2
            if screen(side,sub) == 1
                if side ==1
                subi = subi+1;
                end
                for cond = 1:2
                    [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
                    [~,~,nrepOFF,~] = data_fileguide(R.subname{sub},cond-1);
                    for nr = 1:nrepOFF
                        id = idbank(nr,side,cond);
                        frq = frqbank(band,nr,side,cond); %% THIS USES THE 3rd Band (HB)
                        R.PA.stn_lb_frq = stn_lb_frqbank(nr,side,cond);
                        %                     load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond}  '_' R.siden{side} '_' R.ipsicon])
                        load([R.datapathr R.subname{sub} '\ftdata\virtual_sources_' num2str(nr) '_ROI_' R.condname{cond} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}])
                        
                        %                 load([R.datapathr R.subname{sub} '\ftdata\r' R.subname{sub} '_LCMV_source_' R.condname{cond} 'nrep_' num2str(nr)])
                        % Create data structure
                        Xdata.fsample = R.pp.cont.full.fs;
                        Xdata.label = {vchansave(id).label{[1 refidSave{side,sub}(cond)+1]}};
                        Xdata.trial{1} = vchansave(id).trial{1}([1 refidSave{side,sub}(cond)+1],:);
                        Xdata.time{1} = vchansave(id).time; %{1}
                        
                        databank{side,cond,subi} = Xdata;
                    end
                end
            end
        end
    end
end
a = 1;