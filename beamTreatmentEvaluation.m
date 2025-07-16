 
function beamTreatmentEvaluation()
    % Initialize parameters
    tumorRadius = 500;       % Tumor radius (cm)
    beamEnergy = 160;      % Beam energy (MeV)
    particleType = 'electron'; % 'proton' or 'electron'
    numParticles = 1e12;   % Number of particles
    
    % Generate simulated data
    [doseTumor, doseHealthy] = simulateEnergyDeposition(tumorRadius, beamEnergy, particleType, numParticles);
    
    % Calculate biological effects
    [TCP, NTCP] = calculateBioEffects(doseTumor, doseHealthy, particleType);
    
    % Generate evaluation repor
    printEvaluationReport(TCP, NTCP, particleType, beamEnergy);
end

function [doseTumor, doseHealthy] = simulateEnergyDeposition(radius, energy, type, n)
    % Simulate energy deposition distribution
    depth = 0:0.1:30; % Depth axis (cm)
    
    switch lower(type)
        case 'proton'
            % Proton Bragg peak model
            peakPos = energy/100; % Simplified relationship
            sigma = 0.1*peakPos;
            doseProfile = n * exp(-(depth-peakPos).^2/(2*sigma^2));
            
        case 'electron'
            % Electron exponential attenuation model
            attenuationLength = energy/50; 
            doseProfile = n * exp(-depth/attenuationLength);
            
        otherwise
            error('Unsupported particle type');
    end
    
    % Convert to dose (Gy)
    massTumor = (4/3)*pi*radius^3 * 1.04; % Tumor mass (kg), density≈1.04g/cm³
    doseTumor = trapz(depth(depth<=radius), doseProfile(depth<=radius)) * 1.6e-13 / massTumor;
    
    % Healthy tissue dose (5cm beyond tumor)
    healthyVolume = (4/3)*pi*(radius+5)^3 - (4/3)*pi*radius^3;
    massHealthy = healthyVolume * 1.04;
    doseHealthy = trapz(depth(depth>radius & depth<=radius+5), doseProfile(depth>radius & depth<=radius+5)) * 1.6e-13 / massHealthy;
end

function [TCP, NTCP] = calculateBioEffects(doseTumor, doseHealthy, type)
    % Biological parameters
    alpha_t = 0.3;    % Tumor α (Gy⁻¹)
    beta_t = 0.03;    % Tumor β (Gy⁻²)
    TD50 = 60;        % Healthy tissue TD50 (Gy)
    gamma = 3;        % NTCP slope parameter
    
    % Proton RBE correction
    if strcmpi(type, 'proton')
        RBE = 1.2;
        doseTumor = doseTumor * RBE;
    else
        RBE = 1;
    end
    
    % Tumor Control Probability (LQ model)
    SF = exp(-alpha_t*doseTumor - beta_t*doseTumor^2);
    TCP = exp(-1e9 * SF); % Assuming 1e9 initial cells
    
    % Normal Tissue Complication Probability (Lyman model)
    NTCP = 1 / (1 + (TD50/doseHealthy)^gamma);
end

function printEvaluationReport(TCP, NTCP, type, energy)
    fprintf('\n===== RADIOTHERAPY EVALUATION REPORT =====\n');
    fprintf('Beam Type: %s (%d MeV)\n', upper(type), energy);
    fprintf('Tumor Control Probability (TCP): %.2f%%\n', TCP*100);
    fprintf('Normal Tissue Complication Probability (NTCP): %.2f%%\n', NTCP*100);
    
    if TCP > 90 && NTCP < 5
        fprintf('\nCONCLUSION: Excellent Treatment Plan ✔️\n');
        fprintf('Recommendation: Ready for clinical use\n');
    elseif TCP > 70 && NTCP < 10
        fprintf('\nCONCLUSION: Acceptable Treatment Plan ✔️\n');
        fprintf('Recommendation: Recommended with optimization\n'); 
    else
        fprintf('\nWARNING: Plan Requires Improvement ⚠️\n');
        fprintf('Recommendation: Adjust beam parameters\n');
    end
    fprintf('==========================================\n');
end

