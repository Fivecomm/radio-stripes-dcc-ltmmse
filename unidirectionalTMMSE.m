function T = unidirectionalTMMSE(H, D, Pi, T, l, S_bar, K, L, N, P)
% Recursive routine for unidirectional TMMSE precoding, involving an 
% information matrix V_bar sequentially updated and forwarded from TX 1 to 
% TX L, and statistical information given by Pi. This code has been adapted
% from Python to MATLAB [1].
%
% INPUT:
%   H       = Matrix of N_sim x K x L*N containing the channel estimations
%   D       = DCC matrix for cell-free setup with dimension L x K where 
%             (l,k) is one if AP l serves UE k and zero otherwise
%   Pi      = Matrix containing estimation of statistical parameters
%   T       = List of precoders
%   l       = iterator over L APs
%   S_bar   =
%   K       = Number of User Equipments (UEs)
%   L       = Number of Access Points (APs)
%   N       = Number of antennas per AP
%   P       = Transmitted power for each AP
%
% OUTPUT:
%   T       = List of precoder
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

% Identify the UEs served by AP l
served_UEs = find(D(l, :) == 1);  

if ~isempty(served_UEs)  % Only process if AP l serves some UEs
    Hl = H(:, (l-1)*N+1:l*N);  % Extract channel matrix for AP l

    % Compute precoder using only served UEs
    Fl = pinv(Hl' * Hl + eye(N) / P) * Hl';
    Pl = Hl * Fl;
    Vl = pinv(eye(K) - Pi{l} * Pl) * (eye(K) - Pi{l});
    
    % Ensure T{l} is initialized properly
    if isempty(T{l})
        T{l} = zeros(N, K);  % Preallocate full size
    end

    % Store the computed precoder only for served UEs
    T{l}(:, served_UEs) = Fl(:, served_UEs) * ...
        Vl(served_UEs, served_UEs) * S_bar(served_UEs, served_UEs);
    
    % Recursive call: Only proceed if the next AP serves at least one 
    % overlapping UE
    if l < L  
        next_served_UEs = find(D(l+1, :) == 1);
        overlapping_UEs = intersect(served_UEs, next_served_UEs);

        if ~isempty(overlapping_UEs)  
            % Compute Vl_bar only for the overlapping UEs
            Vl_bar = eye(length(overlapping_UEs)) - ...
                Pl(overlapping_UEs, overlapping_UEs) * ...
                Vl(overlapping_UEs, overlapping_UEs);
            
            % Prepare S_bar_new to forward to the next AP
            S_bar(overlapping_UEs, overlapping_UEs) = ...
                Vl_bar * S_bar(overlapping_UEs, overlapping_UEs);

        else
            % If no overlap, use identity
            S_bar = eye(K);
        end

        % Recursive call with updated S_bar
        T = unidirectionalTMMSE(H, D, Pi, T, l+1, S_bar, K, L, N, P);
        
    end
end

end

