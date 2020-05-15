%%
% 2D Excitation with Spiral trajectory 
% Written by Ali Aghaeifar, 2020, <ali.aghaeifar.mri[at]gmail.com>>
%
%% Timing & Resolution
T   = 2.5e-3;   % Sec, dur of pulse
dt  = 5e-6;     % RF/grad dwell time
t   = 0:dt:T-dt;% time points

fov   = 20;     % cm
res   = 1;      % cm, resolution
kmax  = 1/2/res;% cycles/m, max spatial frequency
N     = ceil(kmax*fov); % number of spiral turns -> kmax / 1 / fov where 1/fov is k-space resolution
gamma = 4257;   % Hz/G

% constant angular rate spiral-in:
k = kmax * (1-t./T) .* exp(1i*2*pi*N*(1-t./T)); % first (1-t./T) indicates spiral-in, second (1-t./T) indicates clock-wise rotation
g = k2g(k, dt); % kspace to gradient conversion

subplot(2,3,1);
plot(k); xlabel('Kx'); ylabel('Ky');
subplot(2,3,2); 
plot(t, real(g), t, imag(g)); ylabel('Gradient, G/cm')

%% build system matrix (2x-oversampled)
xx = -fov/2 : res/2 : fov/2-res/2;
yy = -fov/2 : res/2 : fov/2-res/2;
zz = 0;
[x,y,z] = ndgrid(xx,yy,zz); 
A    = dt*2*pi*gamma*exp( 1i*2*pi*( x(:)*real(k(:)') + y(:)*imag(k(:)') ) );

% build desired pattern
mdes = (2*x.^2 + y.^2 < 5^2);
mdes = IFFT2D( (hamming(size(x,1))*hamming(size(x,1))').^2 .* FFT2D(mdes) );
subplot(2,3,4);
imagesc(abs(mdes)); caxis([-1 1]); title('Desired Pattern');

% design RF
lambda = 1;
b1 = (A'*A + lambda*eye(length(k))) \ (A'*mdes(:)); % = inv(A'*A + lambda*eye(length(k))) * (A'*mdes(:));

subplot(2,3,3);
plot(t, real(b1), t, imag(b1)); ylabel('RF, (G)')

% forward projection
nm = reshape (A*b1, size(x));
subplot(2,3,5);
imagesc(abs(nm)); caxis([-1 1]); title('Estimated Pattern');

%% Simulate it
gr = [real(g(:)), imag(g(:)), zeros(size(g(:)))];
freq = zeros(size(x(:)));
sens = ones(size(b1,2), length(x(:))); 
T1 = 5;
T2 = 2;
[mxcc, mycc, mzcc] = blochCim(b1, gr, dt, T1, T2, freq, [x(:),y(:),z(:)], 0, sens);
mcb = mxcc+1i*mycc;
mcb = reshape(mcb, size(x));

subplot(2,3,6);
imagesc(abs(mcb)); caxis([-1 1]); title('Bloch Simulation');
