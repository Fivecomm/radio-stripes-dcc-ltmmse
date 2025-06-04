function T = UTMMSE(H, D, Pi, K, L, N, P)
% Unidirectional TMMSE precoding, wrapper function for recursive
% computation. This code has been adapted from Python to MATLAB [1].
%
% INPUT:
%   H       = Matrix of N_sim x K x L*N containing the channel estimations
%   D       = DCC matrix for cell-free setup with dimension L x K where 
%             (l,k) is one if AP l serves UE k and zero otherwise
%   Pi      = Matrix containing estimation of statistical parameters
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

% List of precoders
T = cell(1, L);
for l = 1:L
    T{l} = zeros(N, K);
end

% Compute precoders recursively  
T = unidirectionalTMMSE(H,D,Pi,T,1,eye(K),K,L,N,P);

% Convert cell array into matrix
T = cell2mat(T');

end