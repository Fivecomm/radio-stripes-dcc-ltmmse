function C_list = get_parameters_LTMMSE(H_list, D, K, L, N, P, N_sim)
% Monte Carlo estimation of statistical parameters for local TMMSE 
% precoding. This code has been adapted from Python to MATLAB [1].
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

% List of matrices Pi{l} to be computed using statistical information
Pi = cell(1, L);
for l = 1:L
    Pi_l = zeros(K, K); % Initialize Pi_l as a complex matrix
    for n = 1:N_sim
        H = squeeze(H_list(n,:,:)); % Extract the channel matrix
        Hl = H(:, (l-1)*N+1:l*N);   % Extract the channel for AP l

        % Find which UEs AP l serves
        served_UEs = find(D(l,:)==1);

        % Compute per-AP power allocation
        if ~isempty(served_UEs)  % Only compute if AP l serves some UEs
            Hl_served = Hl(served_UEs, :);  % Select assigned UEs

            % Compute statistical parameter using per-UE power
            Pi_l(served_UEs, served_UEs) = ...
                Pi_l(served_UEs, served_UEs) + Hl_served * ...
                pinv(Hl_served' * Hl_served + eye(N) / P) * ...
                Hl_served' / N_sim;
        end

    end
    Pi{l} = Pi_l;
end

% Build linear system of equations using Pi
A = zeros(K*L, K*L) + 1i * zeros(K*L, K*L);
for l = 1:L
    for j = 1:L
        if j == l
            A(K*(l-1)+1:K*l, K*(j-1)+1:K*j) = eye(K);
        else
            A(K*(l-1)+1:K*l, K*(j-1)+1:K*j) = Pi{j};
        end
    end
end

I = zeros(K*L, K) + 1i * zeros(K*L, K);
for l = 1:L
    I(K*(l-1)+1:K*l, :) = eye(K);
end

% Solve system and find optimal coefficients
C = A \ I;

% Convert into list of matrices
C_list = cell(1, L);
for l = 1:L
    C_list{l} = C(K*(l-1)+1:K*l, :);
end

end
