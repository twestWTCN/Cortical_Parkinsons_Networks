clear all; close all
dbstop if error
R.datapathr = 'C:\home\data\TimExtracts310516\';
for sub = 1%:length(R.subnames);
    condNames = {'off','on'};
    for cond = 0:1;
        cd('C:\Users\twest\Documents\Work\PhD\LitvakProject\Coherent_Source_Extract\Scripts')
        ROItable = ROI_Table();
        Ds = dbs_meg_rest_source_cvaraw_TW(R.subnames{sub},cond,[],ROItable{sub}{cond+1},R);
        FTdata = spm2ft_concat(Ds);
        FTdata.history{1} = ['sourceextractloop_ashpointtable_' date];
        save(['C:\Users\twest\Documents\Work\PhD\LitvakProject\Coherent_Source_Extract\Data\' R.subnames{sub} condNames{cond+1}],'FTdata')
    end
end
