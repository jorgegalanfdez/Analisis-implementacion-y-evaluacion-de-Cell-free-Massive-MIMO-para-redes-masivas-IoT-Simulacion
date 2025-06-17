function [SE_MMSE, SE_MR_cent, SINR_MMSE, SINR_MR] = functionComputeSE_uplink(Hhat,H,D,C,tau_c,tau_p,nbrOfRealizations,N,K,p)
%Compute uplink SE for different receive combining schemes using the capacity
%bound in Theorem 5.1 for the centralized schemes and the capacity bound
%in Theorem 5.4 for the distributed schemes. Compute the genie-aided uplink
%SE from Corollary 5.9 for the centralized operation and from Corollary 5.10 
%for the distributed operation. 
%
%INPUT:
%Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel to
%                    UE k in channel realization n.
%H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                    with the true channel realizations. The matrix is
%                    organized in the same way as Hhat.
%D                 = DCC matrix for cell-free setup with dimension L x K 
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%D_small           = DCC matrix for small-cell setup with dimension L x K
%                    where (l,k) is one if AP l serves UE k and zero otherwise
%B                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel estimate 
%                    between AP l and UE k, normalized by noise variance
%C                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix of the channel
%                    estimation error between AP l and UE k,
%                    normalized by noise variance
%tau_c             = Length of coherence block
%tau_p             = Length of pilot sequences
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Number of UEs 
%L                 = Number of APs
%p                 = Uplink transmit power per UE (same for everyone)
%R                 = Matrix with dimension N x N x L x K where (:,:,l,k) is
%                    the spatial correlation matrix between AP l and UE k,
%                    normalized by noise
%pilotIndex        = Vector containing the pilot assigned to each UE
%
%OUTPUT:
%SE_MMSE           = SEs achieved with MMSE combining in (5.11)
%SE_P_MMSE         = SEs achieved with P-MMSE combining in (5.16)
%SE_P_RZF          = SEs achieved with P-RZF combining in (5.18)
%SE_MR_cent        = SEs achieved with centralized MR combining in (5.14)
%SE_opt_L_MMSE     = SEs achieved with opt LSFD in (5.30) and 
%                    L-MMSE combining in (5.29)
%SE_nopt_LP_MMSE   = SEs achieved with n-opt LSFD in (5.41) and 
%                    LP-MMSE combining in (5.39)
%SE_nopt_MR        = SEs achieved with n-opt LSFD in (5.41) and 
%                    local MR combining in (5.32)
%SE_L_MMSE         = SEs achieved with L-MMSE combining in (5.29),
%                    without LSFD   
%SE_LP_MMSE        = SEs achieved with LP-MMSE combining in (5.39),
%                    without LSFD 
%SE_MR_dist        = SEs achieved with local MR combining in (5.32),
%                    without LSFD
%Gen_SE_P_MMSE     = Genie-aided SEs achieved with P-MMSE combining in (5.16)
%Gen_SE_P_RZF      = Genie-aided SEs achieved with P-RZF combining in (5.18)
%Gen_SE_LP_MMSE    = Genie-aided SEs achieved n-opt LSFD in (5.41) and 
%                    LP-MMSE combining in (5.39)
%Gen_SE_MR_dist    = Genie-aided SEs achieved n-opt LSFD in (5.41) and 
%                    local MR combining in (5.32)
%SE_small_MMSE     = SEs achieved with L-MMSE combining for small-cell setup
%Gen_SE_small_MMSE = Genie-aided SEs achieved with L-MMSE combining for 
%                    small-cell setup
%
%
%This Matlab function was developed to generate simulation results to:
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


%Store the N x N identity matrix
eyeN = eye(N);

%Compute the prelog factor assuming only uplink data transmission
prelogFactor = (1-tau_p/tau_c);

%Prepare to store simulation results
SE_MMSE =  zeros(K,1);
SE_MR_cent = zeros(K,1);

SINR_MMSE =  zeros(K,nbrOfRealizations);
SINR_MR = zeros(K,nbrOfRealizations);



%% Go through all channel realizations
for n = 1:nbrOfRealizations
    
    
    %Consider the centralized schemes
        
    %Go through all UEs
    for k = 1:K
                
        %Determine the set of serving APs for UE k
        servingAPs = find(D(:,k)==1); %cell-free setup
        
        %Compute the number of APs that serve UE k
        La = length(servingAPs);
        
        %Determine which UEs that are served by partially the same set of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        %Extract channel realizations and estimation error correlation matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
            
        
        %Compute MMSE combining according to (5.11)
        v = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+eye(La*N))\Hhatallj_active(:,k));
        
        %Compute numerator and denominator of instantaneous SINR in (5.5) for MMSE combining
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        %Update the SE by computing the instantaneous SE for one channel realization according to (5.4)        
        SE_MMSE(k) = SE_MMSE(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        SINR_MMSE(k,n)=numerator/denominator;
        

        
        %Compute centralized MR combining according to (5.14)
        v = Hhatallj_active(:,k);
        
        %Compute numerator and denominator of instantaneous SINR in (5.5)     
        numerator = p*abs(v'*Hhatallj_active(:,k))^2;
        denominator = p*norm(v'*Hhatallj_active)^2 + v'*(p*C_tot_blk+eye(La*N))*v - numerator;
        
        %Update the SE by computing the instantaneous SE for one
        %channel realization according to (5.4)         
        SE_MR_cent(k) = SE_MR_cent(k) + prelogFactor*real(log2(1+numerator/denominator))/nbrOfRealizations;
        SINR_MR(k,n)=numerator/denominator;
        
        
    end
    
end


