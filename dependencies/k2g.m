function [gradients, slewrates] = k2g(kspace, dt, calc_mode)
%
% Convert k-space trajectory to gradients waveform (and slew-rate)
%
% Inputs:
%   kspace    -- kspace trajectory (cycle / cm). n*1, n*2, n*3 real, or n*1 complex.
%   dt        -- time interval between samples (second)
%   calc_mode -- 0 (default): start from kspace = [0,0,0]; 1: start from kspace(1)
% Outputs:
%   gradients -- G/cm
%   slewrates -- G/cm/second
%
% Written by Ali Aghaeifar, 2020, <ali.aghaeifar.mri[at]gmail.com>
% 

if nargin == 0 || isempty(kspace)
    error('Not enough input arguments.')
elseif nargin == 1 || isempty(dt)
    dt = 1e-5; % Gradients raster time in Siemens scanners
    calc_mode = 0;
elseif nargin == 2
    calc_mode = 0;
end
if ~ismatrix(kspace)
    error('k-space must be a vector or matrix');
end
sz = size(kspace);
% To avoid confusion between k-space dimension and samples 
if max(sz) < 4
    error('Minimum samples for k-space is 4.');
end
% Rows are kspace samples
if sz(2)>sz(1) 
    kspace = transpose(kspace); % nonconjugate transpose (Hint: kspace' is conjucate transpose and kspace.' is nonconjucate transpose)
    sz = size(kspace);
end

gamma = 4257;  % Hz/G

if calc_mode == 0
    gradients = -diff( vertcat(zeros(1,sz(2)), kspace) ) / gamma / dt; % also add zero at the beginning to get same length as kspace after diff
elseif calc_mode == 1
    t = transpose(0:dt:(sz(1)-1)*dt); % sz(1) == number of k-space sample points
    g = -diff(kspace) / gamma / dt; % G/cm
    gradients = interp1(t, vertcat(zeros(1,sz(2)), g), t+dt/2, 'linear' ,0); % resample g to get same length as kspace 
else
    error('Third argument must be 0 or 1')
end
slewrates = vertcat(zeros(1,sz(2)), diff(gradients)/dt);

