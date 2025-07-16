%% Electron Scattering Simulation in Gelatin
% Simulates angular distribution of electrons after passing through gelatin layers

clear; clc; close all;

%% Parameters
PARTICLE_TYPE = 'electron'; % Only for electrons
INITIAL_ENERGY = 150; % MeV
N_PARTICLES = 10000; % Number of simulated particles
THICKNESSES = 5:2:35; % Different gelatin thicknesses (mm)
PHANTOM.DENSITY = 1.04; % g/cm³ (gelatin)
PHANTOM.STEP_SIZE = 0.1; % Step size (cm)
PHANTOM.Z = 7.5; % Effective atomic number
PHANTOM.A = 14.0; % Effective atomic mass
PHANTOM.ION_POT = 75e-6; % Mean ionization potential (MeV)

CONST.m_e = 0.510998; % MeV/c² (electron mass)
CONST.c = 299792458e2; % Speed of light (cm/s)

%% Simulation
num_subplots = length(THICKNESSES);
n_cols = ceil(sqrt(num_subplots));
n_rows = ceil(num_subplots / n_cols);
figure('Name', 'Electron Scattering Distribution', 'Position', [100, 100, 1200, 800]);

for idx = 1:num_subplots
    depth = THICKNESSES(idx) / 10; % Convert mm to cm
    
    % Generate initial particle data
    particles.energy = ones(N_PARTICLES, 1) * INITIAL_ENERGY;
    
    % Compute scattering using Molière theory with randomness
    theta_rms = electron_scattering(particles.energy, depth);
    theta_distribution = theta_rms .* randn(N_PARTICLES, 1); % Gaussian spread
    theta_distribution = rad2deg(theta_distribution); % Convert to degrees
    
    % Remove outliers for better visualization
    theta_distribution = theta_distribution(abs(theta_distribution) < 50);
    
    % Reduce number of surviving particles with depth
    survival_fraction = exp(-depth / 2); % Exponential decay approximation (electrons scatter more than protons)
    num_surviving_particles = round(N_PARTICLES * survival_fraction);
    theta_distribution = theta_distribution(1:num_surviving_particles);
    
    % Plot histogram with proper normalization
    subplot(n_rows, n_cols, idx);
    bin_edges = linspace(min(theta_distribution), max(theta_distribution), 50);
    bin_width = bin_edges(2) - bin_edges(1);
    counts = histcounts(theta_distribution, bin_edges, 'Normalization', 'pdf');
    bar(bin_edges(1:end-1) + bin_width/2, counts / sum(counts * bin_width), 'FaceColor', [0.8, 0.2, 0.6], 'EdgeColor', 'none');
    xlabel('Scattering Angle (degrees)');
    ylabel('Probability Density');
    title(sprintf('Thickness = %d mm', THICKNESSES(idx)));
    xlim([-0.5,0.5]);
    ylim([0,8]);
    grid on;
end

%% Function: Electron Scattering (Molière Theory Approximation)
function theta_rms = electron_scattering(E, depth)
    E = max(E, 0.1); % Ensure energy is positive
    theta_rms = (13.6 ./ E) .* sqrt(depth) .* (1 + 0.038*log(depth)); % Molière theory
    theta_rms = deg2rad(theta_rms); % Convert to radians
end
