function dxdt=NRF2p53plusH2O2_DDE(t,x,Z)
%General script to randomly draw a H2O2 generation rate from a lognormal 
%dist with mean = 0.104 and CV = 25%, each simulation runs for 1000 minutes
%This script is called by Call_NRF2p53_TCGA100sims.m

%Random H2O2 generation rate
M=.104;
CV=0.25;

V=(CV*M)^2;

invar= lognrnd(muFunc(M,V),sigmaFunc(M,V));

    function mu = muFunc(M,V)
        mu = log((M^2)/sqrt(V+M^2));
        return
    end


    function sigma = sigmaFunc(M,V)
        sigma = sqrt(log(V/(M^2)+1));
        return
    end



%First order reactions (min^-1), Second order (Cs-1*min^-1), unless
%otherwise noted. 

k_keap1ox = 0.0084; %Oxidation rate of KEAP1 with H2O2 - https://www.ncbi.nlm.nih.gov/pubmed/24634836, https://www.ncbi.nlm.nih.gov/pubmed/20061377
k_keap1red = 6; %Reduction rate of KEAP1 - based on https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2935341/
k_nkon_ETGE = 0.1; %On rate for Nrf2 to Keap1 (based on Khalil et al 2015) https://www.ncbi.nlm.nih.gov/pubmed/25449014
k_nkon_DLG = 0.1; %On rate for Nrf2 to Keap1 (based on Khalil et al 2015) https://www.ncbi.nlm.nih.gov/pubmed/25449014
k_nkD_ETGE = 5E-3; %Dissociation constant for Nrf2 to Keap1 to ETGE motif (ITC data) http://mcb.asm.org/content/27/21/7511.full
k_nkD_DLG = 1; %Dissociation constant for Nrf2 to Keap1 to DLG motif (ITC data) http://mcb.asm.org/content/27/21/7511.full
k_nrf2deg = 0.02; %Nrf2 degradation rate (experimentally determined as detailed in Parameter summary) 
k_nrf2degox = 0.00603; %Nrf2 degradation when keap1 is oxidized - http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030024
k_nrf2cyt = 120; %Rate of Nrf2 moving into nucleus (um^3/min) - based on http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030024
k_nrf2nu = 60; %Rate of Nrf2 moving to cytoplasm (um^3/min) - based on http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.0030024
Ksyn_nrf2 = 0.005; %Synthesis rate of NRF2

%MAF-related kinetic parameters from - https://www.ncbi.nlm.nih.gov/pubmed/16716189
k_nmon = 1; %On rate of Nrf2 to Maf (SPR data) 
k_nmD = 10; %Dissociation constant for Nrf2 to Maf (SPR data)
k_mmon = 1; %On rate of Maf to Maf (SPR data)
k_mmD = 20; %Dissociation constant of Maf to Maf (SPR data)
k_nmAREon = 1.2; %On rate of Nrf2:Maf to ARE sites (SPR data)
k_nmARED = 20; %Dissociation constant of Nrf2:Maf to ARE sites (SPR data)

%Modeling the production of NRF2-target genes (Response, "R") - 
F = .5; %Production rate of ARE proteins 
K_mnARE = 0.2; %Half maximal activation of ARE genes
n = 4; %Hill coefficient 

k_Rred = 0.025; %Reduction of oxidized R to R

Vc = 1000; %Volume of cytosol (um^3)
Vn = 100; %Volume of nucleus (um^3)

%p53-related parameters, from "Recurrent initiation: a mechanism for 
%triggering p53 pulses in response to DNA damage." PMID:18471974 
% Table S2: Model species and parameters
alpha_MPi=5/60; %Mdm2-induced degradation of inactive p53, 1/min
alpha_Pi=2/60;  %inactive p53 degradation, 1/min
alpha_MPa=1.4/60; %Mdm2-induced degradation of active p53, 1/min
alpha_SM=0.5/60; %Signal-dependent inactivation rate of Mdm2, 1/min
alpha_M=1/60; %Mdm2 degradation, 1/min
alpha_i=0.7/60; %Wip1 degradation 1/min
alpha_iS=.01/60; %Wip-induced inactivation of signal, 1/min
beta_P=0.9/60; %p53 production, 1/min
beta_SP=10/60; %signal-induced p53 activation, 1/min
beta_M=0.9/60; %p53-induced Mdm2 production, 1/min
beta_Mi=0.2/60; %p53-independent Mdm2 production, 1/min
beta_i=0.25/60; %p53-induced Wip1 production, 1/min
n_i=4; %Hill coefficient for signal inhibition
n_S=4;  %Hill coefficient of active p53 production by Signal
T_i=0.2; %half-maximal signal inhibition threshold
T_S=8; %half-maximal p53 activation threshold

beta_R_p53=beta_M; %production of R (antioxidants) by p53, modeled as p53 produces MDM2 

%p21-related parameters, from "Fluctuations in p53 Signaling Allow Escape 
%from Cell-Cycle Arrest." PMID:30057196
alpha_mrna_p21=0.22; %Maximum p21 transcription rate, units: 1/min
km_mrna_p21=473.61; %Km of p53-dependent p21 transcription, arbitrary units
beta_mrna_p21=0.0066; %basal p21 mRNA degradation, units 1/min

alpha_prot_p21=1; %p21 protein production rate, units 1/min
beta_prot_p21=0.0023; %Basal p21 protein degradation rate, units 1/min

%ATM-related parameters
k_phos_ATM=.0014; %Phosphorylation rate of ATM in response to hydrogen peroxide, fit to experimental timecourse data (Fig. s9)
k_pATMred=1/60; %basal rate of dephosphorylation (other phosphotases besides Wip1 acting on pATM)
k_syn_ATM=0.003; %Synthesis rate of ATM 
k_deg_ATM=0.0005; %Degradation rate of ATM 

basal_syn_R=0.0006; %Basal R synthesis rate
k_deg_R=0.001; %turnover rate of R reacting with other oxidized proteins in cell not taken into account in this model

k_prod_H2O2=.0052; %Units: uM/min, increased rate from DCIS model      


p53DD=1; %p53 is always mutated in the TNBC models 

if p53DD==1
    beta_M=0; %R19
    beta_i=0; %R23
    alpha_mrna_p21=0; %R28
    beta_R_p53=0; %R30
end

%Define state conditions    
Keap1 = x(1); %Keap1
Keap1_ox = x(2); %oxidized Keap1
Nrf2 = x(3); %Nrf2
Nrf2_Keap1 = x(4); %Nrf2-Keap1
Nrf2_Keap1_ox = x(5); %oxidized NRF2-Keap1
R = x(6); %R (antioxidants)
Nrf2_nucl = x(7); %nuclear NRF2
Maf = x(8); %Maf
Nrf2_Maf = x(9); %Nrf2-Maf heterodimer
Maf_Maf = x(10); %Maf-Maf dimer
Nrf2_Maf_ARE = x(11); %NRF2-Maf bound to antioxidant response element
ARE = x(12); %antioxidant response element
R_ox = x(13); %oxidized R

Pi=x(14); %inactive p53
Pa=x(15); %active p53
M=x(16); %Mdm2
inh=x(17); %inhibitor = Wip1
ATM=x(18); %ATM
p21_mrna=x(19); %p21_mrna
pATM=x(20); %pATM
p21_prot=x(21); %p21_prot

H2O2_in=x(22); %intracellular H2O2

Pa_lag1=Z(:,1);
Pa_lag2=Z(:,2);

%Reaction rates
R1 = k_keap1ox * Keap1 * H2O2_in; %Oxidation of Keap1
R2 = k_keap1red * Keap1_ox * R; %Reduction based on generated reducing proteins
R3 = k_nkon_ETGE * (Nrf2*Keap1 - k_nkD_ETGE*Nrf2_Keap1); %Nrf2 and Keap1 binding through ETGE motif
R3_1 = k_nkon_DLG * (Nrf2*Keap1 - k_nkD_DLG*Nrf2_Keap1); %NRF2 and KEAP1 binding through DLG motif
R4 = Ksyn_nrf2; %synthesis rate of Nrf2
R5 = k_keap1ox * Nrf2_Keap1 * H2O2_in; %Oxidation of Keap1 when in complex with Nrf2
R6 = k_keap1red * Nrf2_Keap1_ox * R; %Reduction of Keap1 when in complex with Nrf2
R7 = k_nrf2deg * Nrf2_Keap1; %Degradation rate of Nrf2
R7_1 = k_nrf2degox * Nrf2_Keap1_ox; %Degradation rate of NRF2 when KEAP1 is oxidized
R8 = k_nrf2cyt * Nrf2 - k_nrf2nu * Nrf2_nucl; %Transport of Nrf2 into nucleus
R9 = k_nmon * (Nrf2_nucl * Maf - k_nmD * Nrf2_Maf); %Nrf2 binding to Maf
R10 = k_mmon * (Maf * Maf - k_mmD * Maf_Maf); %Maf homodimerizing 
R11 = k_nmAREon * (Nrf2_Maf*ARE - k_nmARED * Nrf2_Maf_ARE); %Nrf2:Maf binding to ARE sites
R12 = F * (Nrf2_Maf_ARE^n / (K_mnARE^n + Nrf2_Maf_ARE^n)); %Production of Reducing enzyme "R"
R13 = k_Rred*R_ox; %Reduction of the oxidized reducing enzyme


R14 = beta_P; %p53 production
R15 = alpha_MPi*M*Pi;  %Mdm2-induced degradation of inactive p53
R16 = beta_SP*Pi*((pATM^n_S)/(pATM^n_S+T_S^n_S)); %signal-induced p53 activation
R17 = alpha_Pi*Pi; %inactive p53 degradation

R18 = alpha_MPa*M*Pa; %Mdm2-induced degradation of active p53
R19 = beta_M*Pa_lag1(2); %p53-induced Mdm2 production
R20 = beta_Mi; %p53-independent Mdm2 production
R21 = alpha_SM*pATM*M; %Signal-dependent inactivation of Mdm2
R22 = alpha_M*M; %Mdm2 degradation
R23 = beta_i*Pa_lag2(2); %p53-induced Wip1 production
R24 = alpha_i*inh; %Wip1 degradation
R25 = H2O2_in*k_phos_ATM*ATM;  %H2O2-induced phosphorylation of ATM
R26 = alpha_iS*((inh^n_i)/(inh^n_i+T_i^n_i))*pATM; %Wip1-induced inactivation of signal

R28 = alpha_mrna_p21*(Pa_lag1(2)/(km_mrna_p21+Pa_lag1(2))); %p21 mRNA production
R29 = beta_mrna_p21*p21_mrna; %p21 mRNA degradation
R30 = beta_R_p53*Pa_lag1(2); %Production of reducing enzymes (antioxidants) by p53
R31 = k_pATMred*pATM*R; %Reduction of disulfide bonds by R formed b/w pATM dimer in response to H2O2

R33 = k_syn_ATM; %synthesis rate of ATM
R34 = k_deg_ATM*ATM; %degradation of ATM

R35 = alpha_prot_p21*p21_mrna; %p21 protein production
R36 = beta_prot_p21*p21_prot; %basal p21 protein degradation

R37=basal_syn_R; %basal synthesis rate of antioxidants
R38=k_deg_R*R; %basal oxidation rate of antioxidants (reacting with other oxidized species in the cell besides oxidized KEAP1 and pATM)

R41=k_prod_H2O2; %basal production rate of H2O2
R42=invar*(t<120); %Logical for step-increase in H2O2 generation for 2 hours

%Differential equations
dxdt=zeros(22,1);
dxdt(1) = R2 - R1 - R3 - R3_1; %Keap1
dxdt(2) = R1 - R2; %Keap1 oxidized
dxdt(3) = R4 - R3 - R3_1 - R7 - R7_1 - R8/Vc; %Nrf2 (cytoplasm)
dxdt(4) = R3 + R3_1 - R5 + R6; %Nrf2:Keap1
dxdt(5) = R5 - R6; %Nrf2:Keap1 oxidized
dxdt(6) = R12 - R6 - R2 + R13 -R31+R37-R38+R30; %R
dxdt(7) = R8/Vn - R9; %Nrf2 nucleus
dxdt(8) = -R9 - R10; %Maf
dxdt(9) = R9 - R11; %Maf:Nrf2
dxdt(10) = R10; %Maf:Maf
dxdt(11) = R11; %Nrf2:Maf:ARE
dxdt(12) = -R11; %ARE
dxdt(13)= R6 + R2 - R13+R31; %Oxidized R

dxdt(14) = R14 - R15 - R16 - R17; %inactive p53
dxdt(15) = R16 - R18; %active p53
dxdt(16)= R19 + R20 - R21 - R22; %Mdm2
dxdt(17) = R23 - R24; %Wip1
dxdt(18)= R31-R25+R33-R34; %ATM (unphosphorylated)
dxdt(19) = R28 - R29; %p21_mrna
dxdt(20)= R25-R31-R26; %pATM
dxdt(21) = R35-R36; %p21_prot

dxdt(22)= -R1 - R5-R25+R41+R42;   %intracellular H2O2

end
