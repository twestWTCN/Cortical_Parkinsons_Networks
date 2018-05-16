close all
dbstop if error
addpath('C:\shared')
% R.datapathr = 'C:\home\data\TimExtracts190516\';
% Sub 2 condition ON is broken?

R.subname = {'JP','KB','LM','WB'};
for sub = 1:length(R.subname)
    condNames = {'off','on'};
    for cond = 0:1
        for band = [1 3]
            load(['ROItable_' R.bandname{band}])
            cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
            %             [files, seq, root, details] = dbs_subjects(initials, drug);
            %             [R.origpath initials '\' initials '_R_1_' druglbl{drug+1}]
            if band == 1
                dbs_meg_rest_extract_continuous_headmodel(R.subname{sub}, cond)
            end
%             for side = 1:2
%                 cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
%                 Ds = dbs_meg_rest_source_cvaraw_TW(R.subname{sub},cond,[],ROItable{sub,cond+1,side},R,band);
%                 cd('C:\Users\twest\Documents\Work\GitHub\Cortical_Parkinsons_Networks')
%                 VC_new = spm2ft_concat(Ds);
%                 VC_new.history{1} = ['sourceextractloop_' date];
%                 save([R.datapathr R.subname{sub} '\ftdata\VC_new_ROI_' R.condname{cond+1} '_' R.siden{side} '_' R.ipsicon  '_' R.bandname{band}],'VC_new')
%                 disp([sub cond band side])
%                 clear VC_new
%             end
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