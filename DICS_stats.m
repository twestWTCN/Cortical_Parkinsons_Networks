    tab = []
    rowname = {};
for sub = 1:numel(subname)

load([datapathr subname{sub} '\stats\DICS_ONvOFF'],'dics_cohstattable')
tab = [tab; [dics_cohstattable.ispi.cohstat; dics_cohstattable.contra.cohstat]];
rowname = [rowname; [subname{sub} '_ipsi']; [subname{sub} '_contra']]
cohmax(:,sub) = dics_cohstattable.contra.cohmax;
end

T = table(tab(:,1),tab(:,2),tab(:,3),tab(:,4),tab(:,5),tab(:,6),tab(:,7))
T.Properties.VariableNames = {'ON' 'ON_SEM' 'OFF' 'OFF_SEM' 'ON_OFF' 't_stat' 'P'};
T.Properties.RowNames = rowname;
