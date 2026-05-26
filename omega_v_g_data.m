function data = omega_v_g_data(airspeed, eigSource, numPhysicalModes)
%OMEGA_V_G_DATA Track physical modal branches and compute omega-V-g data.

airspeed = airspeed(:).';

if nargin < 3 || isempty(numPhysicalModes)
    if iscell(eigSource)
        numPhysicalModes = floor(size(eigSource{1}.Aae, 1) / 2);
    else
        numPhysicalModes = floor(size(eigSource, 1) / 2);
    end
end

if iscell(eigSource)
    if numel(eigSource) ~= numel(airspeed)
        error('The models cell array must have one entry for each airspeed value.');
    end
    lambdaBranches = track_physical_mode_branches_from_models(eigSource, numPhysicalModes);
else
    if size(eigSource, 2) ~= numel(airspeed)
        error('eigTraj must have one column for each airspeed value.');
    end
    lambdaBranches = track_physical_mode_branches(eigSource, numPhysicalModes);
end

omega = abs(imag(lambdaBranches));
frequencyHz = omega / (2*pi);

dampingG = nan(size(lambdaBranches));
omegaTol = omega_tolerance(lambdaBranches);
validOmega = omega > omegaTol;
% Stable oscillatory roots have negative g with this convention.
dampingG(validOmega) = 2*real(lambdaBranches(validOmega)) ./ omega(validOmega);

data = struct();
data.airspeed = airspeed;
data.lambda = lambdaBranches;
data.frequencyHz = frequencyHz;
data.dampingG = dampingG;
data.flutter = find_first_flutter_crossing(airspeed, frequencyHz, dampingG);

end

% =========================================================================
function lambdaBranches = track_physical_mode_branches_from_models(models, numPhysicalModes)

nV = numel(models);
lambdaBranches = nan(numPhysicalModes, nV);

firstCandidates = physical_roots_from_model(models{1}, numPhysicalModes);
lambdaBranches(:, 1) = firstCandidates;

for iV = 2:nV
    candidates = physical_roots_from_model(models{iV}, numPhysicalModes);
    lambdaBranches(:, iV) = match_next_branches(lambdaBranches(:, iV-1), candidates);
end

end

% =========================================================================
function roots = physical_roots_from_model(model, numPhysicalModes)

if ~isfield(model, 'Aae') || isempty(model.Aae)
    error('omega-V-g requires models with non-empty Aae matrices.');
end

[eigVectors, eigMatrix] = eig(model.Aae);
eigValues = diag(eigMatrix);

imagTol = omega_tolerance(eigValues);
candidateIdx = find(imag(eigValues) > imagTol);
if numel(candidateIdx) < numPhysicalModes
    error(['Could not find %d positive-frequency roots in Aae at V = %.3f m/s. ', ...
           'Only %d were available.'], numPhysicalModes, model.V, numel(candidateIdx));
end

structuralRows = 1:min(2*numPhysicalModes, size(model.Aae, 1));
candidateVectors = eigVectors(:, candidateIdx);
totalNorm = sum(abs(candidateVectors).^2, 1);
structuralNorm = sum(abs(candidateVectors(structuralRows, :)).^2, 1);
participation = structuralNorm ./ max(totalNorm, eps);

candidateValues = eigValues(candidateIdx);
[~, participationOrder] = sort(participation, 'descend');
roots = candidateValues(participationOrder(1:numPhysicalModes));

[~, frequencyOrder] = sort(abs(imag(roots)), 'ascend');
roots = roots(frequencyOrder);
roots = roots(:);

end

% =========================================================================
function lambdaBranches = track_physical_mode_branches(eigTraj, numPhysicalModes)

nV = size(eigTraj, 2);
lambdaBranches = nan(numPhysicalModes, nV);

firstCandidates = positive_frequency_roots(eigTraj(:, 1));
if numel(firstCandidates) < numPhysicalModes
    error(['Could not find %d positive-frequency roots at the first airspeed. ', ...
           'Only %d were available.'], numPhysicalModes, numel(firstCandidates));
end

lambdaBranches(:, 1) = firstCandidates(1:numPhysicalModes);

for iV = 2:nV
    candidates = positive_frequency_roots(eigTraj(:, iV));
    if numel(candidates) < numPhysicalModes
        error(['Could not find %d positive-frequency roots at airspeed index %d. ', ...
               'Only %d were available.'], numPhysicalModes, iV, numel(candidates));
    end

    lambdaBranches(:, iV) = match_next_branches(lambdaBranches(:, iV-1), candidates);
end

end

% =========================================================================
function roots = positive_frequency_roots(eigValues)

imagTol = omega_tolerance(eigValues);
roots = eigValues(imag(eigValues) > imagTol);
[~, order] = sort(abs(imag(roots)), 'ascend');
roots = roots(order);
roots = roots(:);

end

% =========================================================================
function nextBranches = match_next_branches(previousBranches, candidates)

numBranches = numel(previousBranches);
nextBranches = nan(numBranches, 1);
used = false(numel(candidates), 1);

[~, branchOrder] = sort(abs(imag(previousBranches)), 'ascend');

for k = 1:numBranches
    branchIdx = branchOrder(k);
    availableIdx = find(~used);
    previous = previousBranches(branchIdx);

    scores = branch_distance(candidates(availableIdx), previous);
    [~, bestLocalIdx] = min(scores);
    bestIdx = availableIdx(bestLocalIdx);

    nextBranches(branchIdx) = candidates(bestIdx);
    used(bestIdx) = true;
end

end

% =========================================================================
function scores = branch_distance(candidates, previous)

lambdaScale = max(abs(previous), 1);
omegaScale = max(abs(imag(previous)), 1);

lambdaDistance = abs(candidates - previous) / lambdaScale;
omegaDistance = abs(abs(imag(candidates)) - abs(imag(previous))) / omegaScale;

scores = lambdaDistance + 0.25*omegaDistance;

end

% =========================================================================
function flutter = find_first_flutter_crossing(airspeed, frequencyHz, dampingG)

flutter = empty_flutter_result();
bestV = Inf;

for modeIdx = 1:size(dampingG, 1)
    gBranch = dampingG(modeIdx, :);
    fBranch = frequencyHz(modeIdx, :);

    for iV = 1:numel(airspeed)-1
        g1 = gBranch(iV);
        g2 = gBranch(iV+1);
        f1 = fBranch(iV);
        f2 = fBranch(iV+1);

        if any(isnan([g1, g2, f1, f2]))
            continue
        end

        if g1 <= 0 && g2 >= 0 && (g1 < 0 || g2 > 0)
            V1 = airspeed(iV);
            V2 = airspeed(iV+1);

            if g1 == 0
                Vf = V1;
                ff = f1;
            elseif g2 == 0
                Vf = V2;
                ff = f2;
            else
                Vf = V1 - g1*(V2 - V1)/(g2 - g1);
                ff = f1 + (Vf - V1)*(f2 - f1)/(V2 - V1);
            end
        else
            continue
        end

        if Vf < bestV
            bestV = Vf;
            flutter.hasFlutter = true;
            flutter.Vf = Vf;
            flutter.ff = ff;
            flutter.modeIndex = modeIdx;
            flutter.crossingIndex = iV;
        end
    end
end

end

% =========================================================================
function flutter = empty_flutter_result()

flutter = struct();
flutter.hasFlutter = false;
flutter.Vf = NaN;
flutter.ff = NaN;
flutter.modeIndex = NaN;
flutter.crossingIndex = NaN;

end

% =========================================================================
function tol = omega_tolerance(values)

scale = max(abs(values(:)));
if isempty(scale) || ~isfinite(scale) || scale == 0
    scale = 1;
end

tol = 1e-8 * scale;

end
