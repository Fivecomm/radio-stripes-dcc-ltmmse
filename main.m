% This script simulates the downlink opertation for distributed scenarios
% based on code in [1] and for radio stripes based in [2].%    
%
% REFERENCES:
%   [1] Özlem Tuğfe Demir, Emil Björnson, and Luca Sanguinetti (2021) 
%       “Foundations of User-Centric Cell-Free Massive MIMO”, 
%       Foundations and Trends in Signal Processing: Vol. 14, No. 3-4,
%       pp. 162-472. DOI: 10.1561/2000000109.
%   [2] Lorenzo Miretti, Emil Björnson, David Gesbert, "Team MMSE Precoding
%       with Applications to Cell-free Massive MIMO,” IEEE Transactions on 
%       Wireless Communications, vol. 21, no. 8, pp. 6242-6255, 
%       August 2022.
%       GitHub: https://github.com/emilbjornson/team-MMSE
%
% This is version 1.0 (Last edited: 2025-04-29)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.

%% SETUP

savingFilename = 'results/scenario_05.mat';


L = 64;                     % Number of APs
K = 12;                     % Number of UEs
N = 4;                      % Number of antennas per AP
rAPs = 60;                  % 
rUEs_max = 50;              %
nbrOfSetups = 20000;        % Number of setups with random UE locations
nbrOfRealizations = 100;    % Number of channel realizations per setup
tau_c = 200;                % Length of the coherence block
p = 100;                    % Total uplink transmit power per UE (mW)
rho_tot = 200;              % Total downlink transmit power per AP (mW)
% Decorrelation distance of the shadow fading in (5.43)
decorr = 9;                 
antennaSpacing = 1/2;       % Antenna spacing (in number of wavelengths)
B = 20e6;                   % Communication bandwidth (Hz)
noiseFigure = 7;            % Noise figure (in dB)
% Angular standard deviation in the local scattering model (in radians)
ASD_varphi_deg = 15;        % Azimuth angle
ASD_theta_deg = 15;         % Elevation angle
% Number of pilots per coherence block
tau_p = 10;                 % Not orthogonal in this case
% Compute the prelog factor assuming only downlink data transmission
preLogFactor = (tau_c-tau_p)/tau_c;
% Angular standard deviation in the local scattering model (in radians)
ASD_varphi = deg2rad(ASD_varphi_deg);   %azimuth angle
ASD_theta = deg2rad(ASD_theta_deg);     %elevation angle
% Urban propagation
W = 20;
h = 20;
h_AP = 20;
h_AT = 1.5;
fc = 2;
sigma_sf = 6;


%% INITIALIZATION

SE_LTMMSE_all = zeros(K,nbrOfSetups);   % Local TMMSE (All)
SE_LTMMSE_DCC = zeros(K,nbrOfSetups);   % Local TMMSE (DCC)
SE_UTMMSE_all = zeros(K,nbrOfSetups);   % Unidirectional TMMSE (All)
SE_UTMMSE_DCC = zeros(K,nbrOfSetups);   % Unidirectional TMMSE (DCC)
SE_MMSE_all = zeros(K,nbrOfSetups);     % Cell-free centralized (All)
SE_MMSE_DCC = zeros(K,nbrOfSetups);     % Cell-free centralized (DCC)
SE_LPMMSE_DCC = zeros(K,nbrOfSetups);   % Cell-free distributed (DCC)

R_LTMMSE_all = zeros(K,nbrOfSetups);   % Local TMMSE (All)
R_UTMMSE_all = zeros(K,nbrOfSetups);   % Unidirectional TMMSE (All)
R_MMSE_all = zeros(K,nbrOfSetups);     % Cell-free centralized (All)

time_LTMMSE_DCC = zeros(1, nbrOfSetups);
time_LTMMSE_all = zeros(1, nbrOfSetups);
time_UTMMSE_DCC = zeros(1, nbrOfSetups);
time_UTMMSE_all = zeros(1, nbrOfSetups);
time_MMSE_DCC = zeros(1, nbrOfSetups);
time_LPMMSE_DCC = zeros(1, nbrOfSetups);

%% MAIN LOOP

for n = 1:nbrOfSetups
    % Display simulation progress
    disp(['Setup ' num2str(n) ' out of ' num2str(nbrOfSetups)]);

    % Generate one setup with UEs at random locations
    [gainOverNoisedB,R,pilotIndex,D,~] = generateSetup(L,K,N,rAPs, ...
        rUEs_max,tau_p,ASD_varphi,ASD_theta,B,noiseFigure,sigma_sf, ...
        decorr,W,h,h_AP,h_AT,fc,antennaSpacing);

    %Define the case when all APs serve all UEs
    D_all = ones(L,K);

    % Generate channel realizations, channel estimates, and estimation
    % error correlation matrices for all UEs to the cell-free APs
    [Hhat,H,~,C] = functionChannelEstimates( ...
        R,nbrOfRealizations,L,K,N,tau_p,pilotIndex,p);


    %% RADIO STRIPES (LTMMSE/UTMMSE)
    
    % Equal power allocation
    P = rho_tot/K;

    H_list = permute(Hhat, [2, 3, 1]);

    % Local TMMSE (All)
    tic;
    C_list = get_parameters_LTMMSE(H_list,D_all,K,L,N,P,nbrOfRealizations);
    [SE_LTMMSE_all(:,n), R_LTMMSE_all(:,n)] = simulate_scenario(H_list, ...
        D_all,K,L,N,nbrOfRealizations,P,'LTMMSE',C_list,preLogFactor);
    time_LTMMSE_all(n) = toc;

    % Local TMMSE (DCC)
    tic;
    C_list = get_parameters_LTMMSE(H_list,D,K,L,N,P,nbrOfRealizations);
    SE_LTMMSE_DCC(:,n) = simulate_scenario(H_list,D,K,L,N, ...
        nbrOfRealizations,P,'LTMMSE',C_list,preLogFactor);
    time_LTMMSE_DCC(n) = toc;

    % Unidirectional TMMSE (All)
    tic;
    Pi = get_parameters_UTMMSE(H_list,D_all,K,L,N,P,nbrOfRealizations);
    [SE_UTMMSE_all(:,n), R_UTMMSE_all(:,n)] = simulate_scenario(H_list, ...
        D_all,K,L,N,nbrOfRealizations,P,'UTMMSE',Pi,preLogFactor);
    time_UTMMSE_all(n) = toc;

    % Unidirectional TMMSE (DCC)
    tic;
    Pi = get_parameters_UTMMSE(H_list,D,K,L,N,P,nbrOfRealizations);
    SE_UTMMSE_DCC(:,n) = simulate_scenario(H_list,D,K,L,N, ...
        nbrOfRealizations,P,'UTMMSE',Pi,preLogFactor);
    time_UTMMSE_DCC(n) = toc;

    
    %% CENTRALIZED CELL-FREE MMSE (All)
    tic;
    [SE_MMSE_all(:,n), R_MMSE_all(:,n)] = simulate_scenario(H_list, ...
        D_all,K,L,N,nbrOfRealizations,P,'MMSE','',preLogFactor);
    time_MMSE_all(n) = toc;


    %% CENTRALIZED CELL-FREE MMSE (DCC)
    tic;
    SE_MMSE_DCC(:,n) = simulate_scenario(H_list,D,K,L,N, ...
        nbrOfRealizations,P,'MMSE','',preLogFactor);
    time_MMSE_DCC(n) = toc;


    %% DISTRIBUTED CELL-FREE (LP-MMSE)
    tic;
    % Full uplink power for the computation of precoding vectors using
    % virtual uplink-downlink duality
    p_full = p*ones(K,1);

    % Obtain the expectations for the computation of the terms in
    % (7.13)-(7.15)
    [~,~,~,~,~,~, ...
        signal_LP_MMSE,signal2_LP_MMSE, scaling_LP_MMSE] = ...
        functionComputeExpectations_mex( ...
        Hhat,H,D,C,nbrOfRealizations,N,K,L,p_full);

    % Prepare arrays to store the vectors \tilde{b}_k in (7.25) and 
    % matrices \tilde{C}_{ki} in (7.26)
    bk = zeros(L,K);
    Ck = zeros(L,L,K,K);

    % Go through all UEs
    for k = 1:K
        % Find the APs that serve UE k
        servingAPs = find(D(:,k)==1);
        % The number of APs that serve UE k
        La = length(servingAPs);
        % Compute the vector in (7.25) for UE k (only the non-zero
        % indices correspondig to serving APs are considered)
        bk(1:La,k) = real(vec(signal_LP_MMSE(k,k,servingAPs)))./ ...
            sqrt(scaling_LP_MMSE(servingAPs,k));

        % Go through all UEs
        for i = 1:K
            % Find the APs that serve UE i
            servingAPs = find(D(:,i)==1);
            % The number of APs that serve UE i
            La = length(servingAPs);
            % Compute the matrices in (7.26) (only the non-zero indices
            % are considered)
            if i==k
                Ck(1:La,1:La,k,k) = bk(1:La,k)*bk(1:La,k)';
            else
                Ck(1:La,1:La,k,i) = ...
                    diag(1./sqrt(scaling_LP_MMSE(servingAPs,i))) ...
                    *(vec(signal_LP_MMSE(k,i,servingAPs)) ...
                    *vec(signal_LP_MMSE(k,i,servingAPs))') ...
                    *diag(1./sqrt(scaling_LP_MMSE(servingAPs,i)));

            end

            for j = 1:La
                Ck(j,j,k,i) = signal2_LP_MMSE(k,i,servingAPs(j))/ ...
                    scaling_LP_MMSE(servingAPs(j),i);
            end
        end

    end

    % Take the real part (in the SINR expression,the imaginary terms cancel
    % each other)
    Ck = real(Ck);

    % Equal power allocation
    rho_dist_equal = (rho_tot/tau_p)*ones(L,K);

    % Compute hte square roots of the power allocation coefficients
    % corresponding to (7.24)
    tilrho_equal = sqrt(rho_dist_equal);

    % Go through all UEs
    for k = 1:K
        % Find APs that serve UE k
        servingAPs = find(D(:,k)==1);
        % The number of APs that serve UE k
        La = length(servingAPs);

        % Compute the numerator and denominator of (7.23) for equal
        % and FPA schemes with two different exponents
        numm_equal = abs(bk(1:La,k)'*tilrho_equal(servingAPs,k))^2;
        denomm_equal = 1-numm_equal;

        for i = 1:K
            servingAPs = find(D(:,i)==1);
            La = length(servingAPs);
            denomm_equal = denomm_equal+tilrho_equal(servingAPs,i)'* ...
                Ck(1:La,1:La,k,i)*tilrho_equal(servingAPs,i);

        end

        % Compute SEs using SINRs in (7.23) and Corollary 6.3 for equal
        SE_LPMMSE_DCC(k,n) = preLogFactor*log2(1+numm_equal/denomm_equal);

    end

    time_LPMMSE_DCC(n) = toc;

    %% CENTRALIZED CELL-FREE MMSE (All) - Emil Simulator

    % %Compute the power allocation in (6.36) for distributed precoding
    % rho_dist = zeros(L,K);
    % 
    % gainOverNoise = db2pow(gainOverNoisedB);
    % 
    % for l = 1:L
    % 
    %     %Extract which UEs are served by AP l
    %     servedUEs = find(D_all(l,:)==1);
    % 
    %     %Compute denominator in (6.36)
    %     normalizationAPl = sum(sqrt(gainOverNoise(l,servedUEs)));
    % 
    %     for ind = 1:length(servedUEs)
    % 
    %         rho_dist(l,servedUEs(ind)) = rho_tot*sqrt( ...
    %             gainOverNoise(l,servedUEs(ind)))/normalizationAPl;
    % 
    %     end
    % 
    % end
    % 
    % % Compute SEs for the case where all APs serve all UEs (All)
    % SE_MMSE_all(:,n) = computeSE_MMSE(Hhat,H,D_all,C,preLogFactor, ...
    %     nbrOfRealizations,N,K,L,p,gainOverNoisedB,rho_tot);
    

    %% CENTRALIZED CELL-FREE MMSE (DCC) - Emil Simulator

    % % Compute SEs for DCC case
    % SE_MMSE_DCC(:,n) = computeSE_MMSE(Hhat,H,D,C,preLogFactor, ...
    %     nbrOfRealizations,N,K,L,p,gainOverNoisedB,rho_tot);

end

%% SAVE RESULTS
save(savingFilename, ...
    'SE_MMSE_all',   'SE_MMSE_DCC', 'R_MMSE_all', ...
    'SE_LTMMSE_all', 'SE_LTMMSE_DCC', 'R_LTMMSE_all',...
    'SE_UTMMSE_all', 'SE_UTMMSE_DCC', 'R_UTMMSE_all',...
    'SE_LPMMSE_DCC', ...
    'time_MMSE_all', 'time_MMSE_DCC', ...
    'time_LTMMSE_all', 'time_LTMMSE_DCC', ...
    'time_UTMMSE_all', 'time_UTMMSE_DCC', ...
    'time_LPMMSE_DCC', ...
    '-v7.3');