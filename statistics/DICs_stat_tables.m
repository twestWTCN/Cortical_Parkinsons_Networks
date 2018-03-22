for band = 2%:numel(R.bandname)
    T = table()
    for sub = 1:numel(R.subname)
        load([R.datapathr R.subname{sub} '\stats\DICS_ONvOFF.mat']) % '_' R.bandname{band} '.mat'])
        tablerow =  array2table( dics_cohstattable.contra.cohstat,...
            'VariableNames',{'ON_mean','On_SEM','OFF_mean','OFF_SEM','ON_OFF','T','P'});
        T = [T;tablerow];
    end
end