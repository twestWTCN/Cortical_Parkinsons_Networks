labeln = {'HB Ph-Sync','HB Ctx Amp','HB STN Amp','LB STN Amp'};
for n = 1:4
figure
set(gcf,'Position',[466 129 1101 500])
subplot(1,2,1)
a  = squeeze(dfa_amp_save(1,n,1,:,:,:));
A = a(:); A(A==0) = [];
b  = squeeze(dfa_amp_save(1,n,2,:,:,:));
% b  = squeeze(var_dphidt_save(1,2,:,:));
B = b(:); B(B==0) = [];
barplot_N2(A,B,'DFA Exponent',labeln{n}); ylim([0.5 1.1])

subplot(1,2,2)
a  = squeeze(dfa_amp_save(2,n,1,:,:,:));
A = a(:); A(A==0) = [];
b  = squeeze(dfa_amp_save(2,n,2,:,:,:));
% b  = squeeze(var_dphidt_save(1,2,:,:));
B = b(:); B(B==0) = [];
barplot_N2(A,B,'log Evidence',labeln{n}); ylim([-8 5])
end

savefigure_v2([R.datapathr 'results\images\DFA\'],['DFA_ON_OFF_barplots'],[],[],[]);