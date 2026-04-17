%% Earth -> NEO Rendezvous Porkchop (Falcon Heavy Launch) 
% Grid search over departure date and TOF for Earth -> NEO rendezvous.
% Maps Falcon Heavy launch performance to delivered mass at the asteroid.
%
% Dependencies:
%   LambertArc_ATD2026, EphSS_car, ephNEO, kep2car,
%   date2mjd2000, mjd20002date, getAstroConstants

clear; clc; close all;

%% Setup constants and mission parameters
muSun     = getAstroConstants('Sun', 'mu'); % [km^3/s^2]
secPerDay = 86400;                          % [s/day]
Isp       = 345;                            % [s]
g0        = 9.80665;                        % [m/s^2]
ve        = Isp * g0;                       % exhaust velocity [m/s]
mDryReq   = 2,169;                           % delivered mass requirement at NEO [kg]

%% Target NEO selection
neoFile   = 'All_NEOS_ATA&TD_2018_2019.csv';
targetStr = '2006 SU49';
neoTbl    = readtable(neoFile);

% Find the asteroid in the database
mask = contains(lower(string(neoTbl.full_name)), lower(targetStr));
idx  = find(mask, 1, 'first');

if isempty(idx)
    error('NEO "%s" not found in %s.', targetStr, neoFile);
end

neoName = strtrim(neoTbl.full_name{idx});
neoID   = 11 + idx;

fprintf('\n================ NEO SELECTION ==================\n');
fprintf('Target:    %s\n', targetStr);
fprintf('Matched:   %s (CSV row %d, ephNEO ID %d)\n', neoName, idx, neoID);
fprintf('=================================================\n\n');

%% Falcon Heavy performance (expendable configuration)
C3_vector      = [0   8.7   15    20    25.3  30    40    50    60    80    100];
WetMass_vector = [26500 16800 13600 12200 11300 10400 8800  7400  6200  4200  2700];

max_C3   = 100;
max_vinf = sqrt(max_C3);

%% Define the extended search grid
Dep_Start = date2mjd2000([2028 12 01 0 0 0]); % Start slightly earlier
Dep_End   = date2mjd2000([2029 10 01 0 0 0]);

StepDays      = 2; % finer resolution for the grid
Dep_Grid      = Dep_Start : StepDays : Dep_End;
TOF_Grid_days = 30 : StepDays : 450;              

nD = numel(Dep_Grid);
nT = numel(TOF_Grid_days);

%% Initialize storage matrices for both branches
% Short-way (tm = +1)
FinalMass_SW = NaN(nD, nT);
DVonboard_SW = NaN(nD, nT);
vInfE_SW     = NaN(nD, nT);
vInfN_SW     = NaN(nD, nT);

% Long-way (tm = -1)
FinalMass_LW = NaN(nD, nT);
DVonboard_LW = NaN(nD, nT);
vInfE_LW     = NaN(nD, nT);
vInfN_LW     = NaN(nD, nT);

% Track the best branch (selected by departure V_inf, which drives C3)
FinalMassMatrix = NaN(nD, nT);
DVonboardMatrix = NaN(nD, nT);
vInfE_best      = NaN(nD, nT);
vInfN_best      = NaN(nD, nT);
tm_best         = NaN(nD, nT);

%% Execute grid scan
fprintf('Scanning %d x %d = %d grid points...\n', nD, nT, nD*nT);

for iD = 1:nD
    t_dep = Dep_Grid(iD);
    [r1, vE] = EphSS_car(3, t_dep);
    r1 = r1(:);  
    vE = vE(:);
    
    for iT = 1:nT
        tof_days = TOF_Grid_days(iT);
        t_arr    = t_dep + tof_days;
        tof_sec  = tof_days * secPerDay;
        
        % Get NEO state at arrival
        [r2, vNEO] = neoState(t_arr, neoID, muSun);
        if isempty(r2), continue; end
        
        % Evaluate both Lambert branches independently
        for tm = [+1, -1]
            try
                sol = LambertArc_ATD2026(r1, r2, tof_sec, tm, muSun);
                
                % Skip invalid solutions
                if any(isnan(sol.v1)) || sol.err_km > 1e-2
                    continue
                end
                
                dv_dep = norm(sol.v1 - vE);    % departure V_inf
                dv_arr = norm(vNEO - sol.v2);  % arrival V_inf
                
                % Calculate launcher mass and any excess DV needed
                [mWet, dv_extra] = launcherMass(dv_dep, C3_vector, ...
                    WetMass_vector, max_vinf, max_C3, ve);
                
                if ~isfinite(mWet) || mWet <= 0, continue; end
                
                dv_onboard = dv_extra + dv_arr;
                mFinal = mWet * exp(-(dv_arr * 1000) / ve);
                
                if ~isfinite(mFinal) || mFinal <= 0, continue; end
                
                % Store the results based on the branch
                if tm == +1
                    vInfE_SW(iD,iT)     = dv_dep;
                    vInfN_SW(iD,iT)     = dv_arr;
                    FinalMass_SW(iD,iT) = mFinal;
                    DVonboard_SW(iD,iT) = dv_onboard;
                else
                    vInfE_LW(iD,iT)     = dv_dep;
                    vInfN_LW(iD,iT)     = dv_arr;
                    FinalMass_LW(iD,iT) = mFinal;
                    DVonboard_LW(iD,iT) = dv_onboard;
                end
                
                % Update best solution based on lowest departure V_inf
                if dv_dep < vInfE_best(iD,iT) || ~isfinite(vInfE_best(iD,iT))
                    vInfE_best(iD,iT)      = dv_dep;
                    vInfN_best(iD,iT)      = dv_arr;
                    FinalMassMatrix(iD,iT) = mFinal;
                    DVonboardMatrix(iD,iT) = dv_onboard;
                    tm_best(iD,iT)         = tm;
                end
            catch
                % Silently catch Lambert solver failures
            end
        end
    end
    
    % Simple progress indicator
    if mod(iD, 20) == 0
        fprintf('  %d / %d departure dates done\n', iD, nD);
    end
end
fprintf('Grid scan complete.\n\n');

%% Check optimums and mission feasibility
nGrid   = numel(FinalMassMatrix);
nFinite = sum(isfinite(FinalMassMatrix(:)));
bestPt  = emptyPick();

if nFinite > 0
    [bestPt.mf, idx_best] = max(FinalMassMatrix(:), [], 'omitnan');
    [i, j] = ind2sub(size(FinalMassMatrix), idx_best);
    
    bestPt.dep   = Dep_Grid(i);
    bestPt.tof   = TOF_Grid_days(j);
    bestPt.arr   = bestPt.dep + bestPt.tof;
    bestPt.vInfE = vInfE_best(i,j);
    bestPt.vInfN = vInfN_best(i,j);
    bestPt.tm    = tm_best(i,j);
    bestPt.dv_ob = DVonboardMatrix(i,j);
end

reqMask = isfinite(FinalMassMatrix) & (FinalMassMatrix >= mDryReq);
nFeas   = sum(reqMask(:));
[depEarliest, depLatest, tofMin, tofMax] = feasibleWindow(reqMask, Dep_Grid, TOF_Grid_days);

%% Print top trajectories
printTopTrajectories(FinalMassMatrix, vInfE_best, vInfN_best, tm_best, ...
    DVonboardMatrix, Dep_Grid, TOF_Grid_days, mDryReq, 20, 'BEST BRANCH');

% Print top long-way trajectories
printTopTrajectories(FinalMass_LW, vInfE_LW, vInfN_LW, ...
    -ones(size(FinalMass_LW)), DVonboard_LW, ...
    Dep_Grid, TOF_Grid_days, mDryReq, 20, 'LONG-WAY ONLY');

% Print top short-way trajectories
printTopTrajectories(FinalMass_SW, vInfE_SW, vInfN_SW, ...
    ones(size(FinalMass_SW)), DVonboard_SW, ...
    Dep_Grid, TOF_Grid_days, mDryReq, 20, 'SHORT-WAY ONLY');

%% Generate Plots
X      = Dep_Grid(:).';
Y      = TOF_Grid_days(:).';
Z_mass = FinalMassMatrix.';
Z_dv   = DVonboardMatrix.';

% Best-branch delivered mass (full & zoomed)
plotPorkchop(X, Y, Z_mass, 'Delivered mass at NEO [kg]', ...
    sprintf('Earth-NEO Rendezvous (%s): Delivered Mass — Full', neoName), ...
    mDryReq, bestPt, 'full');

plotPorkchop(X, Y, Z_mass, 'Delivered mass at NEO [kg]', ...
    sprintf('Earth-NEO Rendezvous (%s): Delivered Mass — Zoomed', neoName), ...
    mDryReq, bestPt, 'zoom');

% Best-branch onboard DV (full & zoomed)
plotPorkchop(X, Y, Z_dv, '\Delta v_{onboard} [km/s]', ...
    sprintf('Earth-NEO Rendezvous (%s): Onboard \\Delta v — Full', neoName), ...
    NaN, bestPt, 'full');

plotPorkchop(X, Y, Z_dv, '\Delta v_{onboard} [km/s]', ...
    sprintf('Earth-NEO Rendezvous (%s): Onboard \\Delta v — Zoomed', neoName), ...
    NaN, bestPt, 'zoom', Z_mass, mDryReq);

% Long-way only porkchop
Z_mass_LW = FinalMass_LW.';
plotPorkchop(X, Y, Z_mass_LW, 'Delivered mass at NEO [kg]', ...
    sprintf('Earth-NEO (%s): Delivered Mass — Long-Way Only', neoName), ...
    mDryReq, bestPt, 'full');

% Short-way only porkchop
Z_mass_SW = FinalMass_SW.';
plotPorkchop(X, Y, Z_mass_SW, 'Delivered mass at NEO [kg]', ...
    sprintf('Earth-NEO (%s): Delivered Mass — Short-Way Only', neoName), ...
    mDryReq, bestPt, 'full');

% Branch selection map
figure('Color','w'); hold on; grid on; box on;
branchMap = tm_best.';
branchMap(~isfinite(branchMap)) = 0;

imagesc(X, Y, branchMap);
% Color mapping: grey=none, blue=short, red=long
colormap([0.8 0.8 0.8; 0.2 0.4 0.8; 0.8 0.2 0.2]);  
cb = colorbar;
cb.Label.String = 'Branch (blue=short, red=long)';
xlabel('Departure date'); 
ylabel('Time of Flight [days]');
title(sprintf('Lambert Branch Selection Map (%s)', neoName));

ax = gca; 
formatDateAxis(ax, ax.XLim);
hold off;

%% Console Summary
fprintf('\n================== SUMMARY ==================\n');
fprintf('Target NEO:        %s (ephNEO ID %d)\n', neoName, neoID);
fprintf('Grid points:       %d\n', nGrid);
fprintf('Finite solutions:  %d (%.1f%%)\n', nFinite, 100*nFinite/nGrid);

nSW = sum(isfinite(FinalMass_SW(:)));
nLW = sum(isfinite(FinalMass_LW(:)));
fprintf('Short-way solutions: %d\n', nSW);
fprintf('Long-way solutions:  %d\n', nLW);

if isfinite(bestPt.mf)
    fprintf('\nOptimum (max delivered mass, selected by min dep V_inf):\n');
    fprintf('  mf     = %.1f kg\n', bestPt.mf);
    fprintf('  Dep    = %s | TOF %.1f d | Arr %s\n', ...
        localDateStr(bestPt.dep), bestPt.tof, localDateStr(bestPt.arr));
    fprintf('  v_inf  = %.3f km/s (Earth) | %.3f km/s (NEO) | C3 = %.2f km^2/s^2\n', ...
        bestPt.vInfE, bestPt.vInfN, bestPt.vInfE^2);
    fprintf('  DV_ob  = %.3f km/s\n', bestPt.dv_ob);
    fprintf('  Branch = %s\n', branchStr(bestPt.tm));
else
    fprintf('\nNo finite Lambert solutions found.\n');
end

fprintf('\nRequirement (mf >= %.0f kg):\n', mDryReq);
fprintf('  Feasible points: %d (%.1f%%)\n', nFeas, 100*nFeas/nGrid);

if nFeas > 0
    fprintf('  Dep window:      %s to %s (step %d d)\n', ...
        localDateStr(depEarliest), localDateStr(depLatest), StepDays);
    fprintf('  TOF range:       %.0f to %.0f days\n', tofMin, tofMax);
    
    tmFeas = tm_best(reqMask);
    fprintf('  Branch split:    short = %d, long = %d\n', ...
        sum(tmFeas == +1, 'omitnan'), sum(tmFeas == -1, 'omitnan'));
end
fprintf('=============================================\n');

%% Local Helper Functions
function [r, v] = neoState(t, id, muSun)
    try
        kep = ephNEO(t, id);
        st  = kep2car(kep, muSun);
        r = st(1:3)';  r = r(:);
        v = st(4:6)';  v = v(:);
    catch
        r = [];  v = [];
    end
end

function [mWet, dv_extra] = launcherMass(vInfE, C3_vec, mass_vec, max_vinf, max_C3, ve)
    dv_extra = 0;
    if vInfE <= max_vinf
        mWet = interp1(C3_vec, mass_vec, vInfE^2, 'linear');
    else
        mWet0    = interp1(C3_vec, mass_vec, max_C3, 'linear');
        dv_extra = vInfE - max_vinf;
        mWet     = mWet0 * exp(-(dv_extra * 1000) / ve);
    end
end

function out = emptyPick()
    out = struct('dep',NaN,'tof',NaN,'arr',NaN,'mf',NaN, ...
                 'vInfE',NaN,'vInfN',NaN,'tm',NaN,'dv_ob',NaN);
end

function [depE, depL, tofMin, tofMax] = feasibleWindow(reqMask, Dep_Grid, TOF_Grid_days)
    depE = NaN; depL = NaN; tofMin = NaN; tofMax = NaN;
    if ~any(reqMask(:)), return; end
    
    [iD, iT] = find(reqMask);
    depE   = min(Dep_Grid(iD));
    depL   = max(Dep_Grid(iD));
    tofMin = min(TOF_Grid_days(iT));
    tofMax = max(TOF_Grid_days(iT));
end

function printTopTrajectories(MassMatrix, vInfE_mat, vInfN_mat, tm_mat, DV_mat, ...
                              Dep_Grid, TOF_Grid_days, mDryReq, nTop, label)
    validMask = isfinite(MassMatrix) & (MassMatrix >= mDryReq);
    if ~any(validMask(:))
        fprintf('\n======== TOP TRAJECTORIES (%s): NONE FEASIBLE ========\n', label);
        return
    end
    
    tempMass = MassMatrix(:);
    tempMass(~validMask(:)) = -Inf;
    [~, globalIdx] = sort(tempMass, 'descend');
    nTop = min(nTop, sum(validMask(:)));
    
    fprintf('\n================ TOP %d TRAJECTORIES (%s) ================\n', nTop, label);
    fprintf('%-4s %-12s %-12s %-7s %-11s %-10s %-10s %-10s %-10s\n', ...
        'Rank','Departure','Arrival','TOF','C3(km2/s2)','Vinf_d','Vinf_a','Mass(kg)','Branch');
    fprintf('--------------------------------------------------------------------------------------------\n');
    
    for k = 1:nTop
        gIdx = globalIdx(k);
        [ii, jj] = ind2sub(size(MassMatrix), gIdx);
        
        depMJD = Dep_Grid(ii);
        arrMJD = depMJD + TOF_Grid_days(jj);
        
        depDate = mjd20002date(depMJD);
        arrDate = mjd20002date(arrMJD);
        
        vinf_d = vInfE_mat(ii,jj);
        vinf_a = vInfN_mat(ii,jj);
        c3val  = vinf_d^2;
        tmVal  = tm_mat(ii,jj);
        
        fprintf('%-4d %04d-%02d-%02d  %04d-%02d-%02d  %-7.0f %-11.2f %-10.3f %-10.3f %-10.1f %-10s\n', ...
            k, ...
            depDate(1), depDate(2), depDate(3), ...
            arrDate(1), arrDate(2), arrDate(3), ...
            TOF_Grid_days(jj), ...
            c3val, vinf_d, vinf_a, ...
            MassMatrix(ii,jj), ...
            branchStr(tmVal));
    end
    fprintf('=============================================================================\n');
end

function plotPorkchop(X, Y, Z, cbLabel, titleStr, mReq, bestPt, mode, varargin)
    figure('Color','w'); hold on; grid on; box on;
    contourf(X, Y, Z, 25, 'LineStyle','none');
    
    cb = colorbar;
    cb.Label.String = cbLabel;
    xlabel('Departure date');
    ylabel('Time of Flight [days]');
    title(titleStr);
    
    if isfinite(mReq)
        [C, h] = contour(X, Y, Z, [mReq mReq], 'k', 'LineWidth', 2);
        try
            clabel(C, h, 'FontSize', 10, 'Color','k', ...
                   'FontWeight','bold', 'LabelSpacing', 450);
        catch
        end
    end
    
    if isfinite(bestPt.mf)
        plot(bestPt.dep, bestPt.tof, 'kp', 'MarkerFaceColor','w', ...
             'MarkerSize', 11, 'LineWidth', 1.2);
    end
    
    if strcmp(mode, 'zoom')
        if nargin >= 10
            Z_ref = varargin{1};  
            thresh = max(0.4*varargin{2}, 200);
        else
            Z_ref = Z;  
            thresh = max(0.4*mReq, 200);
        end
        
        zmask = isfinite(Z_ref) & (Z_ref >= thresh);
        if any(zmask(:))
            [iY, iX] = find(zmask);
            xPad = 0.10 * (X(max(iX)) - X(min(iX)));
            yPad = 0.10 * (Y(max(iY)) - Y(min(iY)));
            xlim([X(min(iX)) - xPad, X(max(iX)) + xPad]);
            ylim([Y(min(iY)) - yPad, Y(max(iY)) + yPad]);
        end
    end
    
    ax = gca;
    formatDateAxis(ax, ax.XLim);
    hold off;
end

function formatDateAxis(ax, xlims)
    xt = linspace(xlims(1), xlims(end), 8);
    ax.XTick = xt;
    xtlbl = strings(size(xt));
    for k = 1:numel(xt)
        d = mjd20002date(xt(k));
        xtlbl(k) = sprintf('%04d-%02d-%02d', d(1), d(2), d(3));
    end
    ax.XTickLabel = xtlbl;
    ax.XTickLabelRotation = 25;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.MinorGridAlpha = 0.15;
    ax.GridAlpha      = 0.25;
    ax.TickDir        = 'out';
    ax.LineWidth      = 1;
end

function s = branchStr(tm)
    if ~isfinite(tm),  s = 'N/A';
    elseif tm > 0,     s = 'short-way';
    else,              s = 'long-way';
    end
end

function s = localDateStr(mjd)
    d = mjd20002date(mjd);
    s = sprintf('%04d-%02d-%02d', d(1), d(2), d(3));
end
