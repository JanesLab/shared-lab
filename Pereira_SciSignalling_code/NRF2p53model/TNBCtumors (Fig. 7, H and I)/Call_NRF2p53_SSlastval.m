function [a]=Call_NRF2p53_SSlastval()

% Initial conditions
y0 = zeros(22,1); 

vals=readtable('TCGA_FCs_medians.csv','ReadRowNames',true);

parfor i=1:width(vals)
    y0=[
        1.74985959173128 * vals{'KEAP1',vals.Properties.VariableNames{i}}
        0.000127644027974550 * vals{'KEAP1',vals.Properties.VariableNames{i}}
        0.0717899022306803 * vals{'NFE2L2',vals.Properties.VariableNames{i}}
        0.249994528315562 * vals{'NFE2L2',vals.Properties.VariableNames{i}} *vals{'KEAP1',vals.Properties.VariableNames{i}}
        1.82359251686784e-05 * vals{'NFE2L2',vals.Properties.VariableNames{i}} *vals{'KEAP1',vals.Properties.VariableNames{i}}
        2.03866217137325  * vals{'Rmed',vals.Properties.VariableNames{i}}
        0.143579627519012 * vals{'NFE2L2',vals.Properties.VariableNames{i}}
        3.35052078657992* vals{'Mafmed',vals.Properties.VariableNames{i}}
        0.0481075208107867 * vals{'NFE2L2',vals.Properties.VariableNames{i}}* vals{'Mafmed',vals.Properties.VariableNames{i}}
        0.561299604083659* vals{'Mafmed',vals.Properties.VariableNames{i}}
        0.0400720885256285 * vals{'NFE2L2',vals.Properties.VariableNames{i}}* vals{'Mafmed',vals.Properties.VariableNames{i}}
        16.6599279114744
        0.160074235636796* vals{'Rmed',vals.Properties.VariableNames{i}}
        
        0.303124246906332 * vals{'TP53',vals.Properties.VariableNames{i}}
        0.000204534293877026 * vals{'TP53',vals.Properties.VariableNames{i}}
        0.193782720693684 * vals{'MDM2',vals.Properties.VariableNames{i}}
        4.58048789188940e-05 * vals{'PPM1D',vals.Properties.VariableNames{i}}
        5.96009077727445 * max(vals{'ATM',vals.Properties.VariableNames{i}},vals{'CHEK2',vals.Properties.VariableNames{i}})
        9.02099752811825e-06* vals{'CDKN1A',vals.Properties.VariableNames{i}}
        0.0652598606320947 * max(vals{'ATM',vals.Properties.VariableNames{i}},vals{'CHEK2',vals.Properties.VariableNames{i}}) %pATM
        0.00329239437042024* vals{'CDKN1A',vals.Properties.VariableNames{i}}
        0.106221780406244];

tspan = [0 5000];
lags=[42, 75]; %lags, minutes (p53-induced expression of Mdm2 delay, p53-induced expression of Wip1 delay)
sol = dde23(@NRF2p53noH2O2_DDE,lags,y0, tspan)

sol_lastval{i}=sol.y(:,end)

end 

save 'SS_lastvalues.mat',sol_lastval
