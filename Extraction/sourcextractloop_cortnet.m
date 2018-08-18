function sourcextractloop_cortnet(R)
cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
if nargin<1
    R = makeHeader_SubCort_Cort_Networks();
end
FRESH = 1; %!!!
for sub = 1:length(R.subname)
        if FRESH ==1
            eval(['! rmdir ' R.datapathr R.subname{sub} ' /s /q'])
            mkdir([R.datapathr R.subname{sub}])
        end
    for cond = 0:1
        for breg = 1:length(R.bregname)
            cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
            if breg == 1
                dbs_meg_rest_extract_continuous_headmodel(R.subname{sub}, cond)
            end
            for side = 1:2
                ROItable = [{['ipsi_' R.siden{side}(1) R.ref_list{1}(1:3) '_' R.siden{side}(1) R.bregname{breg}]}...
                    {'all'},{R.bandef(R.bregband{breg},:)},{R.bregROI{breg}(side,:)}];
                cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
                Ds = dbs_meg_rest_source_cvaraw_CPN(R.subname{sub},cond,[],ROItable,R,breg,R.siden{side}(1),R.pp.cont.thin.bp);
                cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
                VC_new = spm2ft(Ds);
                VC_new.history{1} = ['sourceextractloop_' date];
                mkdir([R.datapathr R.subname{sub} '\ftdata\'])
                save([R.datapathr R.subname{sub} '\ftdata\VC_new_ROI_' R.condname{cond+1} '_' R.siden{side} '_' R.ipsicon  '_' R.bregname{breg}],'VC_new')
                disp([sub cond breg side])
                clear VC_new
            end
        end
    end
end

    function FTdata = spm2ft(Ds)
        FTdata.label = Ds.chanlabels;  % cell-array containing strings, Nchan X 1
        
        FTdata.fsample = Ds.fsample; % sampling frequency in Hz, single number
        [m n t] = size(Ds(:,:,:));
        B= squeeze(mat2cell(Ds(:,:,:), m, n, ones(1,t)))';
        FTdata.trial = B; % cell-array containing a data matrix for each trial (1 X Ntrial), each data matrix is    Nchan X Nsamples
        FTdata.time = repmat({repmat(Ds.time,size(Ds,1),1)},1,size(Ds,3));   % cell-array containing a time axis for each trial (1 X Ntrial), each time axis is a 1 X Nsamples vector
        %         FTdata.trialinfo = repmat([subnames{sub} '_Rest1_' condName{cond}],1,1); % this field is optional, but can be used to store trial-specific information, such as condition numbers, reaction times, correct responses etc. The dimensionality is N x M
    end
end