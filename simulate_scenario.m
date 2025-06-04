function [SE,R] = simulate_scenario(H_list, D, K, L, N, N_sim, P, ...
    precoder, parameters, preLogFactor)
% Numerically evaluate the performance of a given precoding scheme, using 
% the MSE lower bound. Since the schemes are optimal, the lower bound 
% coincides with the UatF bound (with unitary UL TX power), or with the 
% hardening bound after applying the UL-DL duality principle. This code has
% been adapted from Python to MATLAB [1].
%
% INPUT:
%   H_list      = Matrix of N_sim x K x L*N with the channel estimations
%   D           = DCC matrix for cell-free setup with dimension L x K where 
%                 (l,k) is one if AP l serves UE k and zero otherwise
%   K           = Number of User Equipments (UEs)
%   L           = Number of Access Points (APs)
%   N           = Number of antennas per AP
%   P           = Transmitted power for each AP
%   N_sim       = Number of channel realizations
%   precoder    = Type of precodig:
%                 - 'MMSE' for centralized MMSE
%                 - 'LTMMSE' for local TMMSE
%                 - 'UTMMSE' for unidirectional TMMSE
%   parameters  = Two possible parameters:
%                 - C_list: optimal coefficients for LTMMSE
%                 - Pi: estimation of statistical parameters for UTMMSE
%   preLogFactor= prelog factor
%
% OUTPUT:
%   SE = Achievable spectral efficiency (bit/s/Hz)
%
%
% REFERENCES:
%   [1] Lorenzo Miretti, Emil Björnson, David Gesbert, “Team MMSE Precoding
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

%% FUNCTION

% Initialize variables
signal_MMSE = zeros(1, K);
interf_MMSE = zeros(1, K);
MSE = zeros(1, K);
sigma2 = 1; % Default noise variance

SE = zeros(1, K); % Initialize the SE array

for n = 1:N_sim
    H = squeeze(H_list(n,:,:));

    % Compute precoders
    switch precoder
        case 'MMSE'
            T = pinv(H' * H + eye(N*L) / P) * H';
        case 'LTMMSE'
            T = LTMMSE(H, D, parameters, K, L, N, P);
        case 'UTMMSE'
            T = UTMMSE(H, D, parameters, K, L, N, P);
    end

    % Update MSE estimate for all users
    for k = 1:K
        e = zeros(K, 1); 
        e(k) = 1; % Identity vector for UE k
        
        % Compute MSE using the given formula
        MSE(k) = MSE(k) + (norm(e - H * T(:, k))^2 + ...
            norm(T(:, k))^2 / P) / N_sim;
    end
    
    % Compute SINR and SE for all users
    for k = 1:K
        % Find APs that serve UE k
        servingAPs = find(D(:,k)==1);

        % Calculate the column indices for the corresponding APs
        colIndices = (servingAPs - 1) * N + (1:N);
        colIndices_k = colIndices(:)';

        % Signal power (only from serving APs)
        signal_power = abs(H(k,colIndices_k) * T(colIndices_k, k))^2;
        
        % Interference power (sum over all users except k)
        interference_power = 0;
        for j = 1:K
            if j ~= k
                % Find APs serving UE j
                servingAPs_j = find(D(:,j)==1);

                % Calculate the column indices for the corresponding APs
                colIndices = (servingAPs_j - 1) * N + (1:N);
                colIndices_j = colIndices(:)';
                
                % Compute interference contribution
                interference_power = interference_power + ...
                    abs(H(k, colIndices_j) * T(colIndices_j, j))^2;
            end
        end
        
        % Noise power
        noise_power = norm(T(colIndices_k, k))^2 / P;
        
        % SINR for user k
        SINR_k = signal_power / (interference_power + noise_power);
        
        % Update SE (averaging over simulations)
        SE(k) = SE(k) + preLogFactor * log2(1 + SINR_k) / N_sim;
    end

end

% Compute achievable rates (SE)
R = -log2(MSE);

end
