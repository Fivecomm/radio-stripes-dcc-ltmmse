function SE = computeSE_MMSE(Hhat,H,D,C,preLogFactor,nbrOfRealizations, ...
    N,K,L,p,gainOverNoisedB,rho_tot)
% Compute downlink SE for the MMSE precoding scheme using the capacity 
% bound in Theorem 6.1 for the centralized schemes and the capacity bound 
% in Corollary 6.3 for the distributed schemes. Compute the genie-aided 
% downlink SE from Corollary 6.6 for the centralized and the distributed
% operations. 
%
% INPUT:
%   Hhat              = Matrix with dimension L*N  x nbrOfRealizations x K
%                       where (:,n,k) is the estimated collective channel 
%                       to UE k in channel realization n.
%   H                 = Matrix with dimension L*N  x nbrOfRealizations x K
%                       with the true channel realizations. The matrix is
%                       organized in the same way as Hhat.
%   D                 = DCC matrix for cell-free setup with dimension L x K 
%                       where (l,k) is one if AP l serves UE k and zero 
%                       otherwise
%   C                 = Matrix with dimension N x N x L x K where (:,:,l,k)
%                       is the spatial correlation matrix of the channel
%                       estimation error between AP l and UE k, normalized 
%                       by noise variance
%   preLogFactor      = Prelog factor
%   nbrOfRealizations = Number of channel realizations
%   N                 = Number of antennas per AP
%   K                 = Number of UEs 
%   p                 = Uplink transmit power per UE (same for everyone)
%   gainOverNoisedB   = Matrix with dimension L x K where (l,k) is the 
%                       channel gain (normalized by the noise variance) 
%                       between AP l and UE k
%   rho_tot           = Maximum allowed transmit power for each AP 
%
% OUTPUT:
%   SE                = SEs achieved with MMSE precoding in (6.16)
%
%
% This Matlab function is a modified version of the
% functionComputeSE_downlink.m function in [1].
%
% REFERENCES:
%   [1] Özlem Tuğfe Demir, Emil Björnson, and Luca Sanguinetti (2021) 
%       “Foundations of User-Centric Cell-Free Massive MIMO”, 
%       Foundations and Trends in Signal Processing: Vol. 14, No. 3-4,
%       pp. 162-472. DOI: 10.1561/2000000109.
%
% This is version 1.0 (Last edited: 2025-04-29)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.
%%

signal_MMSE = zeros(K,1);
interf_MMSE = zeros(K,1);
scaling_MMSE = zeros(K,1);
portionScaling_MMSE = zeros(L,K);

%% Compute scaling factors for precoding

%Go through all channel realizations
for n = 1:nbrOfRealizations

    %Go through all UEs
    for k = 1:K
        
        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        La = length(servingAPs);
        
        %Determine which UEs that are served by partially the same set
        %of APs as UE k, i.e., the set in (5.15)
        servedUEs = sum(D(servingAPs,:),1)>=1;
        
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        C_tot_blk_partial = zeros(N*La,N*La);
        
        for l = 1:La
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape( ...
                Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = ...
                sum(C(:,:,servingAPs(l),:),4);
            C_tot_blk_partial((l-1)*N+1:l*N,(l-1)*N+1:l*N) = ...
                sum(C(:,:,servingAPs(l),servedUEs),4);
        end
        
        %Compute MMSE, P-MMSE, and P-RZF precoding
        V_MMSE = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+ ...
            eye(La*N))\Hhatallj_active(:,k));
        
        %Compute scaling factor by Monte Carlo methods
        scaling_MMSE(k) = scaling_MMSE(k) + sum(abs(V_MMSE).^2,1)/ ...
            nbrOfRealizations;
        
        %Go through all the serving APs
        for l=1:La
            
            %Extract the portions of the centralized precoding vectors 
            V_MMSE2 = V_MMSE((l-1)*N+1:l*N,:);
            
            %Compute realizations of the terms inside the expectations
            %of the signal and interference terms in Theorem 6.1
            portionScaling_MMSE(servingAPs(l),k) = portionScaling_MMSE( ...
                servingAPs(l),k)+sum(abs(V_MMSE2).^2,1)/nbrOfRealizations;
           
        end
    end

end

%Normalize the norm squares of the portions for the normalized centralized 
%precoders
portionScaling_MMSE = portionScaling_MMSE./repmat(scaling_MMSE.',[L 1]);

%% Power allocation

%The parameters for the scalable centralized downlink power allocation in
%(7.43)
upsilon = -0.5;
kappa = 0.5;

%Compute the power allocation coefficients for centralized precoding
%according to (7.43)
rho_MMSE = functionCentralizedPowerAllocation( ...
    K,gainOverNoisedB,D,rho_tot,portionScaling_MMSE,upsilon,kappa);

%% Go through all channel realizations
for n = 1:nbrOfRealizations

    %Matrix to store Monte-Carlo results for this realization
    interf_MMSE_n = zeros(K,K);

    %Go through all UEs
    for k = 1:K

        %Determine the set of serving APs
        servingAPs = find(D(:,k)==1);
        
        La = length(servingAPs);
        
        %Extract channel realizations and estimation error correlation
        %matrices for the APs that involved in the service of UE k
        Hallj_active = zeros(N*La,K);
        Hhatallj_active = zeros(N*La,K);
        C_tot_blk = zeros(N*La,N*La);
        
        for l = 1:La
            Hallj_active((l-1)*N+1:l*N,:) = reshape(...
                H((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            Hhatallj_active((l-1)*N+1:l*N,:) = reshape(...
                Hhat((servingAPs(l)-1)*N+1:servingAPs(l)*N,n,:),[N K]);
            C_tot_blk((l-1)*N+1:l*N,(l-1)*N+1:l*N) = ...
                sum(C(:,:,servingAPs(l),:),4);
        end        

        %Compute MMSE combining
        w = p*((p*(Hhatallj_active*Hhatallj_active')+p*C_tot_blk+ ... 
            eye(La*N))\Hhatallj_active(:,k));
            
        %Apply power allocation
        w = w*sqrt(rho_MMSE(k)/scaling_MMSE(k));
        
        %Compute realizations of the terms inside the expectations
        %of the signal and interference terms in Theorem 6.1
        signal_MMSE(k) = signal_MMSE(k) + (Hallj_active(:,k)'*w)/ ...
            nbrOfRealizations;
        interf_MMSE_n(:,k) = interf_MMSE_n(:,k) + Hallj_active'*w;
    end

    %Compute interference power in one realization
    interf_MMSE = interf_MMSE + sum(abs(interf_MMSE_n).^2,2)/ ...
        nbrOfRealizations;

end

SE = preLogFactor*real(log2(1+(abs(signal_MMSE).^2) ./ ...
    (interf_MMSE - abs(signal_MMSE).^2 + 1)));

end

