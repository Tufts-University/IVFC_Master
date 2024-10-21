%% TH_052924_PressureCalculator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% file:NV_050223_PressureCalculator.m and TH_052924_PressureCalculator.m
% ***Description***:
% This function calculates the needed input pressure for the pressure 
% driven flow based on viscosity, number of channels, and channel gap. 
% Written By: Nilay Vora (nvora01@tufts.edu)
% Date Written: 05/02/2023
% Modifying Author:Taras Hanulia
% Date Modified:05/29/2024
% the presure at the channel shoud be divided to number of channel. Reff
% modification.
% Latest Revision: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialization
clear 
clc
%% Viscosity 
% Water: 0.89 cP
% Blood: 4 cP
% Media: 0.93 cP
% Rat blood 1.53 cP
%% Unit Conversions:
% 1 uL = 1 mm^3
% 1 cP = 1e-3 Pa*s
% 1 min = 60 sec
% 1 Pa = 1e-2 mbar
% 1 um = 1e-3 mm
%% Default Values
r_tube = 0.38/2; % mm (from tubing specifications)
P_atm = 1013.25; % mbar
Q = 2; % uL/min (Desired Flow Rate)
Q = Q./60; % mm^3/sec
height = 32;% um (channel height)
height = height./1000; % mm
L_split = 1.46; % mm (Distance from center of punch to channel split)
L_channels = 10; % mm (Length of Channels)
W_channels = 0.03; % mm (channel width)
%% Equations 
R_circ = @(eta,L) (8.*eta*L)./(pi*r_tube.^4);
R_rect = @(eta,L,width) (28.4*eta*L)./...
                        (height.^4); 
R_rect2 = @(eta,L,width) (12*eta*L)./...
                        (width*height.^3.*(1-0.63*(height/width))); 
width = @(n_channel,w_gap) W_channels*n_channel + (n_channel-1)*w_gap; % um
P0 = @(R_eff) Q.*R_eff;
%% Obtain values
prompt = {'Enter number of channels:', 'Enter channel gap size (um):',...
          'Enter length of inlet tubing (mm):',...
          'Enter length of outlet tubing (mm):',...
          'Enter sample Viscosity (cP):'};
definput = {'1','0','360','90', '0.89'}; % default values
dims = [1 50];
dlgtitle = 'Pressure Calculator';
answer = inputdlg(prompt,dlgtitle,dims,definput);
n_channels = str2double(answer{1}); % Number of channels
w_gap = str2double(answer{2})./1000; % Width of gap (mm)
L_in = str2double(answer{3}); % Tubing length into device (mm)
L_out = str2double(answer{4}); % Tubing length out of device (mm)
eta = str2double(answer{5}) * 1e-3 * 1e-2; % Viscosity (milibar * s)
%% Calculate Pressure
width_final = width(n_channels,w_gap);
if width_final>2*height
    factor_R = @(eta,L_split,width_final) R_rect2(eta,L_split,width_final);
else
    factor_R = @(eta,L_split,width_final) R_rect(eta,L_split,width_final);
end
R_eff = @(n_channel) R_circ(eta,L_in) + R_circ(eta,L_out) + ...
        2.* factor_R(eta,L_split,width_final) + ...
        factor_R(eta,L_channels,W_channels)./n_channel;
R_final = R_eff(n_channels);
P_final = n_channels.*P0(R_final);
disp('-------------------------------------------------------------------')
disp('Pressure needed for specified parameters =')
disp([num2str(P_final), ' milibars'])
disp('-------------------------------------------------------------------')
