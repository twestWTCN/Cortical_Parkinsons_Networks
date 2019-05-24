function [] = group_DICS_imagemean(R)
for sub = 1:numel(R.subname)
    %     Tmri = ft_read_mri(['C:\Users\Tim\Documents\Work\Cortical_Networks\MRI_processing\template_MRI\single_subj_T1_1mm.nii'],'dataformat','nifti_spm');
    %     Tmri = ft_convert_units(Tmri,'cm');
    Tmri = ft_read_mri([R.datapathr R.subname{sub} '\MRI\orig\r' R.subname{sub} '.nii'],'dataformat','nifti_spm');
    Tmri = ft_convert_units(Tmri,'cm');
    for cond = 1:2
        [datafileN,pp_mark,nrep,senscheck] = data_fileguide(R.subname{sub},cond-1);
        for nr = nrep
            parfor refN = 1:numel(R.ref_list)
                ref_chan = R.ref_list{refN};
                a = R.subname{sub}; b = R.condname{cond};
                source = [];
                dataload = load([R.datapathr a '\ftdata\r' a '_DICSv2_source' b 'nrep_' num2str(nr) '_' ref_chan],'source')
                source = dataload.source;
                T = [1  0   0    0
                    0   1   0    0
                    0   0   1    2.2
                    0   0   0    1];
                source = ft_transform_geometry(T,source);
                
                cfg            = [];
                cfg.downsample = 2;
                cfg.parameter = 'avg.coh';
                sourceInt  = ft_sourceinterpolate(cfg, source , Tmri);
                source_avg_dics_cell{refN,cond,sub} = sourceInt;
                source_avg_dics(:,refN,cond,sub) = sourceInt.coh(:);
                source_dims(:,refN,cond,sub) = sourceInt.dim;
                [refN,cond,sub]
            end
        end
    end
end

save([R.datapathr 'results\images\groupDICsresults'],'source_avg_dics','source_dims')
save([R.datapathr 'results\images\groupDICsresults_sourcecell'],'source_avg_dics_cell','-v7.3')


