function pathloss = functionPathlossNLoS(W,h,h_AP,h_AT,fc,d)
% Estimate the path loss baed on the "NLoS" propagation models (urban,
% suburban, rural) specified in [1].
%
% INPUT:
%   W           = Street width (m)
%   h           = Average building height (m)
%   h_AP        = AP antenna height (m)
%   h_AT        = UE AT antenna height (m)
%   fc          = carrier frequency (GHz)
%   d           = distance between transmitter and receiver antenna (m)
%
% OUTPUT:
%   pathloss    = path loss (large-scale fading coefficient) (dB)
%
%
% REFERENCES:
%   [1] M. Series, “Guidelines for evaluation of radio interface 
%       technologies for IMT-Advanced," Report ITU M.2135-1, 2009. 
%
% This is version 1.0 (Last edited: 2023-10-10)
%
% License: This code is licensed under the GPLv2 license. If you in any way
% use this code for research that results in publications, please cite our
% monograph as described above.

pathloss = 161.04 - 7.1 * log10(W) + 7.5 * log10(h) ...
    - (24.37 - 3.7 * (h / h_AP).^2) * log10(h_AP) ...
    + (43.42 - 3.1 * log10(h_AP)) * (log10(d) - 3) ...
    + 20 * log10(fc) - (3.2 * (log10(11.75 * h_AT)).^2 - 4.97);

end