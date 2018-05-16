function grandAverageSpectra_Screened(R)
load([R.datapathr 'results\spectral\screen.mat'])
close all
for band = [1 3]
    for sub = 1:numel(R.subname)
        for side = 1:2
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_power_' R.ipsicon '_'  R.bandname{band}],'powsave')
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_coh_' R.ipsicon '_'  R.bandname{band}],'cohsave','frqsave')
            load([R.datapathr R.subname{sub} '\ftdata\ROI_analy\ROIvoxel_bank_' R.ipsicon '_'  R.bandname{band}],'frqbank','idbank','stn_lb_frqbank','refsave')
            [~,~,nrep(1),~] = data_fileguide(R.subname{sub},0);
            [~,~,nrep(2),~] = data_fileguide(R.subname{sub},1);
            id(1) = idbank(nrep(1),side,1);
            id(2) = idbank(nrep(2),side,2);
            
            
            if screen(side,sub) == 1
                for ch = 1:2
                    ON = powsave{nrep(1),side,1};
                    OFF = powsave{nrep(2),side,2};
                    ON = squeeze(ON(ch,:,:)); if size(ON,1)<size(ON,2); ON = ON'; end
                    OFF = squeeze(OFF(ch,:,:));if size(OFF,1)<size(OFF,2); OFF = OFF'; end
                    frq = frqsave{nrep(1),side,1}(:,1)'; powsubgrand{1,ch,side,sub} = ON(:,id(1));
                    frq = frqsave{nrep(2),side,2}(:,1)'; powsubgrand{2,ch,side,sub} = OFF(:,id(2));
                end
                ON =  cohsave{nrep(1),side,1}; cohsubgrand{1,side,sub} = ON(:,id(1));
                OFF = cohsave{nrep(2),side,2}; cohsubgrand{2,side,sub} = OFF(:,id(2));
                
            end
            
        end
    end
    spectralplots_groups(R.datapathr,powsubgrand,cohsubgrand,frq,R.titular,R.bandname{band})
end