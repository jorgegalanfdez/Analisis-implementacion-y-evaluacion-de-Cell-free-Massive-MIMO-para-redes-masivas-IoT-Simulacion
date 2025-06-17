%This Matlab script can be used to reproduce Figures 5.4(a) and 5.6(a) in the monograph:
%
%Ozlem Tugfe Demir, Emil Bjornson and Luca Sanguinetti (2021),
%"Foundations of User-Centric Cell-Free Massive MIMO", 
%Foundations and Trends in Signal Processing: Vol. 14: No. 3-4,
%pp 162-472. DOI: 10.1561/2000000109
%
%This is version 1.0 (Last edited: 2021-01-31)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%monograph as described above.

%Empty workspace and close figures
close all;
clear;


%% Define simulation setup

%Number of Monte-Carlo setups
nbrOfSetups = 100;

%Number of channel realizations per setup
nbrOfRealizations = 100;

%Number of APs 
L = 100;

%Number of antennas per AP
N = 1;

%Number of UEs in the network
K = 20; % K = [20 40 60 80 100 120 140 160 180 200]

%Length of coherence block
tau_c = 200;

%Length of pilot sequences
tau_p = 20;

%Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(15);  %azimuth angle
ASD_theta = deg2rad(15);   %elevation angle

%% Propagation parameters

%Total uplink transmit power per UE (mW)
p = 100; % p = [200 150 100 90 80 70 60 50 40 30 20 10 1]

%Prepare to save simulation results
SE_MMSE_original = zeros(K,nbrOfSetups); %MMSE (All)
SINR_MMSE_original = zeros(K,nbrOfRealizations,nbrOfSetups);
SE_MR_original = zeros(K,nbrOfSetups); %MR (All)
SINR_MR_original = zeros(K,nbrOfRealizations,nbrOfSetups);


SE_MMSE_DCC = zeros(K,nbrOfSetups); %MMSE (DCC)
SINR_MMSE_DCC = zeros(K,nbrOfRealizations,nbrOfSetups);
SE_MR_DCC = zeros(K,nbrOfSetups); %MR (DCC)
SINR_MR_DCC = zeros(K,nbrOfRealizations,nbrOfSetups);

%% Go through all setups
for n = 1:nbrOfSetups
    
    %Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);
    
    %Generate one setup with UEs and APs at random locations
    [gainOverNoisedB,R,pilotIndex,D,~] = generateSetup(L,K,N,tau_p,1,0,ASD_varphi,ASD_theta);
    
    %Generate channel realizations with estimates and estimation
    %error correlation matrices
    [Hhat,H,~,C] = functionChannelEstimates(R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);
    
    
    
    %% Original Cell-Free Massive MIMO
    
    %Define the case when all APs serve all UEs
    D_all = ones(L,K);
    
    %Compute SE using combiners and results in Section 5 for centralized
    %and distributed uplink operations for the case when all APs serve all UEs
    [SE_MMSE_all, SE_MR_cent_all, SINR_MMSE_all, SINR_MR_all]=functionComputeSE_uplink(Hhat,H,D_all,C,tau_c,tau_p,nbrOfRealizations,N,K,p);
    
    %Save SE values for setups
    SE_MMSE_original(:,n) = SE_MMSE_all;
    SINR_MMSE_original(:,:,n) = SINR_MMSE_all;
    SE_MR_original(:,n) = SE_MR_cent_all;
    SINR_MR_original(:,:,n) = SINR_MR_all;
    


    %% Cell-Free Massive MIMO with DCC
    
    %Compute SE using combiners and results in Section 5 for centralized
    %and distributed uplink operations for DCC
    [SE_MMSE, SE_MR_cent, SINR_MMSE, SINR_MR]=functionComputeSE_uplink(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,p);
    
    %Save SE values for setups
    SE_MMSE_DCC(:,n) =  SE_MMSE;
    SINR_MMSE_DCC(:,:,n) = SINR_MMSE;
    SE_MR_DCC(:,n) =  SE_MR_cent;
    SINR_MR_DCC(:,:,n) = SINR_MR;
    
    %Remove large matrices at the end of analyzing this setup
    clear Hhat H B C R;
    
end

Setup_results_file=['CFmMIMO_MR_MMSE_UL_K_' num2str(K) '_p_' num2str(p) '.mat'];
save(Setup_results_file,'SE_MMSE_original','SINR_MMSE_original','SE_MR_original','SINR_MR_original','SE_MMSE_DCC','SINR_MMSE_DCC','SE_MR_DCC','SINR_MR_DCC');


%% Plot simulation results
% figure;
% hold on; box on;
% set(gca,'fontsize',16);
% 
% plot(sort(SE_MMSE_original(:)),linspace(0,1,K*nbrOfSetups),'k-','LineWidth',2);
% plot(sort(SE_MMSE_DCC(:)),linspace(0,1,K*nbrOfSetups),'r-.','LineWidth',2);
% plot(sort(SE_MR_original(:)),linspace(0,1,K*nbrOfSetups),'b--','LineWidth',2);
% plot(sort(SE_MR_DCC(:)),linspace(0,1,K*nbrOfSetups),'k:','LineWidth',3);
% 
% xlabel('Spectral efficiency [bit/s/Hz]','Interpreter','Latex');
% ylabel('CDF','Interpreter','Latex');
% legend({'MMSE (All)','MMSE (DCC)','P-MMSE (DCC)','P-RZF (DCC)','MR (DCC)'},'Interpreter','Latex','Location','SouthEast');
% xlim([0 12]);

