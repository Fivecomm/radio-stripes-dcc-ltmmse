function Pi = get_parameters_UTMMSE(H_list, D, K, L, N, P, N_sim)
% Monte Carlo estimation of statistical parameters for unidirectional TMMSE
% precoding. Iterative implementation, could have been implemented 
% recursively similarly to the unidirectionalTMMSE routine.
% This code has been adapted from Python to MATLAB [1].
%
% INPUT:
%   H_list  = Matrix of N_sim x K x L*N containing the channel estimations
%   D       = DCC matrix for cell-free setup with dimension L x K where 
%             (l,k) is one if AP l serves UE k and zero otherwise
%   K       = Number of User Equipments (UEs)
%   L       = Number of Access Points (APs)
%   N       = Number of antennas per AP
%   P       = Transmitted power for each AP
%   N_sim   = Number of channel realizations
%
% OUTPUT:
%   C_list = Matrix containing optimal coefficients
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

% List of matrices Pi[l] to be computed using statistical information
Pi = cell(1, L);
for l = 1:L
    Pi{l} = zeros(K, K);
end

% Iterative computation from TX L to TX 1
for l = L:-1:2
    % Estimate auxiliary statistical quantities
    E_PS = zeros(K, K);
    E_Sbar = zeros(K, K);

    for n = 1:N_sim
        H = squeeze(H_list(n,:,:));
        Hl = H(:, (l-1)*N+1:l*N);    % Extract AP l’s channel matrix

        % Find UEs assigned to AP l
        served_UEs = find(D(l, :) == 1);  

        if ~isempty(served_UEs)  % Only process if AP l serves UEs
            Hl_served = Hl(served_UEs, :);  % Select only assigned UEs
            
            % Compute statistical matrices only for served UEs
            Pl = Hl_served * pinv(Hl_served' * Hl_served + eye(N) / P) *...
                Hl_served';
            Vl = pinv(eye(length(served_UEs)) - ...
                Pi{l}(served_UEs, served_UEs) * Pl) * ...
                (eye(length(served_UEs)) - ...
                Pi{l}(served_UEs, served_UEs));
            Vl_bar = pinv(eye(length(served_UEs)) - Pl * ...
                Pi{l}(served_UEs, served_UEs)) * ...
                (eye(length(served_UEs)) - Pl);
            
            % Accumulate expectations for assigned UEs
            E_PS(served_UEs, served_UEs) = ...
                E_PS(served_UEs, served_UEs) + Pl * Vl / N_sim;
            E_Sbar(served_UEs, served_UEs) = ...
                E_Sbar(served_UEs, served_UEs) + Vl_bar / N_sim;
        end

    end
    % Update Pi[l-1] using Pi[l] and auxiliary statistical quantities
    Pi{l-1} = E_PS + Pi{l} * E_Sbar;
end

end