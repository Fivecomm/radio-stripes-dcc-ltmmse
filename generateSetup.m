function [gainOverNoisedB,R,pilotIndex,D,D_small,APpositions, ...
    UEpositions,distances] = generateSetup(L, K, N, rAPs, rUEs_max, ...
    tau_p, ASD_varphi, ASD_theta, B, noiseFigure, sigma_sf, decorr, W, ...
    h, h_AP, h_AT, fc, antennaSpacing)
% This function is an adaptation of the function that generates 
% realizations of the simulation setup described in [1, Section 5.3].
%
% INPUT:
%   L               = Number of APs per setup
%   K               = Number of UEs in the network
%   N               = Number of antennas per AP
%   rAPs            = Radius where APs are equally distributed (m)
%   rUEs_max        = Radius within which UEs are distributed (m)
%   tau_p           = Number of orthogonal pilots
%   ASD_varphi      = Angular standard deviation in the local scattering 
%                     model for the azimuth angle (in radians)
%   ASD_theta       = Angular standard deviation in the local scattering 
%                     model for the elevation angle (in radians)
%   B               = Communication bandwidth (Hz)
%   noiseFigure     = Noise figure (in dB)
%   sigma_sf        = Standard deviation of the shadow fading in (5.43)
%   decorr          = Decorrelation distance of the shadow fading in (5.43)
%   W               = Street width (m)
%   h               = Average building height (m)
%   h_AP            = AP antenna height (m)
%   h_AT            = UE AT antenna height (m)
%   fc              = carrier frequency (GHz)
%   antennaSpacing  = Define the antenna spacing (in number of wavelengths)
%
% OUTPUT:
%   gainOverNoisedB = Matrix with dimension L x K where element (l,k) is 
%                     the channel gain (normalized by the noise variance)
%                     between AP l and UE k in setup n
%   R               = Matrix with dimension N x N x L x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between AP l and UE k in setup n, normalized by noise
%   pilotIndex      = Matrix with dimension K x nbrOfSetups containing the
%                     pilots assigned to the UEs
%   D               = DCC matrix with dimension L x K x nbrOfSetups where 
%                     (l,k,n) is one if AP l serves UE k in setup n and
%                     zero otherwise for cell-free setup
%   D_small         = DCC matrix with dimension L x K x nbrOfSetups where
%                     (l,k,n) is one if AP l serves UE k in setup n and
%                     zero otherwise for small-cell setup
%   APpositions     = Vector of length L with the AP locations, where the
%                     real part is the horizontal position and the 
%                     imaginary part is the vertical position
%   UEpositions     = Vector of length K with UE positions, measured in the
%                     same way as APpositions
%   distances       = Matrix with same dimension as gainOverNoisedB 
%                     containing the distances in meter between APs and UEs
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
% use this code for research that results in publications, please cite [1]
% as described above.


%% Define simulation setup

rng('shuffle')

%Prepare to save results
gainOverNoisedB = zeros(L,K);
R = complex(zeros(N,N,L,K));
distances = zeros(L,K);
pilotIndex = zeros(K,1);
D = zeros(L,K);
D_small = zeros(L,K);

masterAPs = zeros(K,1); %the indices of master AP of each UE k 

%% Simulate setup

%The square length for the wrap around is the same as the APs circle
squareLength = rAPs;

%Compute noise power (in dBm)
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Positions of APs in a circle of radio r_APs equally spaced
theta_APs = linspace(0, 2*pi, L+1)';
APpositions = rAPs .* (cos(theta_APs(1:L)) + sin(theta_APs(1:L)) * 1i);

%Random UEs locations with uniform distribution inside radio r_UEs
rUEs = rand(1, K) * rUEs_max;
theta_UEs =  rand(1, K) * 2 * pi;
UEpositions = rUEs .* (cos(theta_UEs) + sin(theta_UEs) * 1i);

%Compute alternative AP locations by using wrap around
wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
wrapVertical = wrapHorizontal';
wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + ...
    repmat(wrapLocations,[L 1]);

%Prepare to store shadowing correlation matrix
shadowCorrMatrix = sigma_sf^2*ones(K,K);
shadowAPrealizations = zeros(K,L);

%Add UEs
for k = 1:K
    
    %Generate a random UE location in the area
    UEposition = UEpositions(k);
    
    %Compute distances assuming that the APs are 10 m above the UEs
    [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped - ...
        repmat(UEposition,size(APpositionsWrapped))),[],2);
    distanceVertical = h_AP - h_AT;
    distances(:,k) = sqrt(distanceVertical^2+distanceAPstoUE.^2);
    
    %If this is not the first UE
    if k-1>0
        
        %Compute distances from the new prospective UE to all other UEs
        shortestDistances = zeros(k-1,1);
        
        for i = 1:k-1
            shortestDistances(i) = ...
                min(abs(UEposition - UEpositions(i) + wrapLocations));
        end
        
        %Compute conditional mean and standard deviation necessary to
        %obtain the new shadow fading realizations, when the previous
        %UEs' shadow fading realization have already been generated.
        %This computation is based on Theorem 10.2 in "Fundamentals of
        %Statistical Signal Processing: Estimation Theory" by S. Kay
        newcolumn = sigma_sf^2*2.^(-shortestDistances/decorr);
        term1 = newcolumn'/shadowCorrMatrix(1:k-1,1:k-1);
        meanvalues = term1*shadowAPrealizations(1:k-1,:);
        stdvalue = sqrt(sigma_sf^2 - term1*newcolumn);
        
    else %If this is the first UE
        
        %Add the UE and begin to store shadow fading correlation values
        meanvalues = 0;
        stdvalue = sigma_sf;
        newcolumn = [];
        
    end
    
    %Generate the shadow fading realizations
    shadowing = meanvalues + stdvalue*randn(1,L);
    
    %Compute the channel gain divided by noise power including shadowing
    pathloss = functionPathlossNLoS(W,h,h_AP,h_AT,fc,distances(:,k));
    gainOverNoisedB(:,k) = -pathloss + shadowing' - noiseVariancedBm;
    
    %Update shadowing correlation matrix and store realizations
    shadowCorrMatrix(1:k-1,k) = newcolumn;
    shadowCorrMatrix(k,1:k-1) = newcolumn';
    shadowAPrealizations(k,:) = shadowing;
    
    %Store the UE position
    UEpositions(k) = UEposition;  
    
    %Determine the master AP for UE k by looking for AP with best
    %channel condition
    [~,master] = max(gainOverNoisedB(:,k));
    D(master,k) = 1;
    masterAPs(k) = master;
    
    %Assign orthogonal pilots to the first tau_p UEs according to
    %Algorithm 4.1
    if k <= tau_p
        
        pilotIndex(k) = k;
        
    else %Assign pilot for remaining UEs
        
        %Compute received power to the master AP from each pilot
        %according to Algorithm 4.1
        pilotinterference = zeros(tau_p,1);
        
        for t = 1:tau_p
            
            pilotinterference(t) = ...
                sum(db2pow(gainOverNoisedB(master,pilotIndex(1:k-1)==t)));
            
        end
        
        %Find the pilot with the least receiver power according to
        %Algorithm 4.1
        [~,bestpilot] = min(pilotinterference);
        pilotIndex(k) = bestpilot;
        
    end
            
    %Go through all APs
    for l = 1:L
        
        %Compute nominal angle between UE k and AP l
        angletoUE_varphi = ... %azimuth angle
            angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); 
        angletoUE_theta = ... %elevation angle
            asin(distanceVertical/distances(l,k));  
        %Generate spatial correlation matrix using the local
        %scattering model in (2.18) and Gaussian angular distribution
        %by scaling the normalized matrices with the channel gain
        % R(:,:,l,k) = db2pow(gainOverNoisedB(l,k)) * ...
        %     functionRlocalscattering(N,angletoUE_varphi, ...
        %     angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);
        R(:,:,l,k) = db2pow(gainOverNoisedB(l,k)) * ...
            functionRlocalscattering_mex(N,angletoUE_varphi, ...
            angletoUE_theta,ASD_varphi,ASD_theta,antennaSpacing);

    end
    
end

%Each AP serves the UE with the strongest channel condition on each of
%the pilots in the cell-free setup
for l = 1:L
    
    for t = 1:tau_p
        
        pilotUEs = find(t==pilotIndex(:));
        [~,UEindex] = max(gainOverNoisedB(l,pilotUEs));
        D(l,pilotUEs(UEindex)) = 1;
       
    end
    
end

%Determine the AP serving each UE in the small-cell setup according to
%(5.47) by considering only the APs from the set M_k for UE k, i.e.,
%where D(:,k,n) is one.
for k=1:K
    
    tempmat = -inf*ones(L,1);
    tempmat(D(:,k)==1,1) = gainOverNoisedB(D(:,k)==1,k);
    [~,servingAP] = max(tempmat);
    D_small(servingAP,k) = 1;
    
end

