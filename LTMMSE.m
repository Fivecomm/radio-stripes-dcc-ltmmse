function T = LTMMSE(H, D, C_list, K, L, N, P)
% Local TMMSE precoding.
% This code has been adapted from Python to MATLAB [1].
%
% INPUT:
%   H       = Matrix of N_sim x K x L*N containing the channel estimations
%   D       = DCC matrix for cell-free setup with dimension L x K where 
%             (l,k) is one if AP l serves UE k and zero otherwise
%   C_list  = Matrix containing optimal coefficients
%   K       = Number of User Equipments (UEs)
%   L       = Number of Access Points (APs)
%   N       = Number of antennas per AP
%   P       = Transmitted power for each AP
%
% OUTPUT:
%   T       = Precoder
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

% Initialize T
T = zeros(L * N, K);

for l = 1:L
    % Find UEs assigned to AP l
    served_UEs = find(D(l, :) == 1);  % UEs served by AP l

    if ~isempty(served_UEs)  % Only process if AP l serves some UEs
        Hl = H(:, (l-1)*N + 1:l*N);  % Channel matrix for AP l
        % Select only the rows corresponding to the UEs served by AP l
        Hl_served = Hl(served_UEs, :);
        % Compute TMMSE precoding for the UEs served by AP l
        Fl = pinv(Hl_served' * Hl_served + eye(N) / P) *  Hl_served';
        % Update precoder T for the served UEs
        T((l-1)*N + 1:l*N, served_UEs) = ...
            Fl * C_list{l}(served_UEs, served_UEs);
    end
end

end

