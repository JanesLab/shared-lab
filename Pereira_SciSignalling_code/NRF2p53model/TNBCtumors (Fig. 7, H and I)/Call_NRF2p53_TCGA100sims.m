function [a]=Call_NRF2p53_TCGA100sims()

nloop=100;

% Initial conditions
y0 = zeros(22,1);
load('SS_lastvalues.mat','sol_lastval')

%Inputs
tspan = [0 1000];
lags=[42, 75]; %lags, minutes (p53-induced expression of Mdm2 delay, p53-induced expression of Wip1 delay)

parfor i=1:length(sol_lastval)
    y0=sol_lastval{i};
    soln={};
    for jz = 1:nloop
        sol = dde23(@NRF2p53plusH2O2_DDE,lags,y0, tspan)
        
        soln{jz}.allspecies=sol.y;
        soln{jz}.time=sol.x;
    end
    
   m=matfile(sprintf('Tumor_%d.mat',i),'writable',true)
   m.x=soln; 
    
end

end


