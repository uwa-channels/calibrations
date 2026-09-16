function run_matlab()
% RUN_MATLAB Run the MATLAB implementation and write its results to artifacts/.
%
% Reads the same config.json and the same probe.mat as the Python runner, so
% neither language owns a second copy of the parameters or draws its own
% random numbers for the probe.
%
% Produces
%    artifacts/matlab_replay.mat   y, y_same_start, y_resamp
%    artifacts/matlab_noise.mat    w, beta
%
% ``start`` in config.json is the Python (0-based) index, so MATLAB is given
% ``start + 1`` throughout.  A second, deliberately misaligned replay is run at
% the same integer: it is one sample of fs_delay out, and the comparer checks
% that it looks clearly worse.  That is the harness testing itself -- a
% comparison that cannot see a one-sample offset would pass everything.
%
% Author: Zhengnan Li
% Email : uwa-channels@ofdm.link
% License: MIT

here = fileparts(mfilename('fullpath'));
addpath(getenv_or(fullfile(here, '..', 'replay_matlab', 'src'), 'UWA_MATLAB_SRC'));

cfg = jsondecode(fileread(fullfile(here, 'config.json')));
artifacts = fullfile(here, 'artifacts');
if ~exist(artifacts, 'dir'); mkdir(artifacts); end

data = getenv_or(here, 'UWA_CHANNELS_CACHE');
channel_path = fullfile(data, cfg.data.channel_file);
noise_path = fullfile(data, cfg.data.noise_file);

%% Probe: the same samples the Python runner used
probe = load(fullfile(artifacts, 'probe.mat'));
x = double(probe.input(:));

%% Replay
channel = load(channel_path);
fs = cfg.replay.fs;
fs_delay = channel.params.fs_delay;
array_index = cfg.replay.array_index(:).' + 1;   % config is 0-based
start = cfg.replay.start;

[p, q] = rat(fs_delay/fs);
down = resample(x, p, q);
y_resamp = resample(down, q, p);

y = replay(x, fs, array_index, channel, start+1);
y_same_start = replay(x, fs, array_index, channel, start);

% The same replay with the phase trajectory removed, so that h_hat, the spline
% interpolation and the two resamplings are the only things acting.
static = rmfield(channel, 'phi_hat');
y_static = replay(x, fs, array_index, static, start+1);

fprintf('matlab replay: y %dx%d, control %d, p/q = %d/%d\n', ...
    size(y, 1), size(y, 2), numel(y_resamp), p, q);

save(fullfile(artifacts, 'matlab_replay.mat'), 'y', 'y_same_start', ...
    'y_static', 'y_resamp', 'start', 'array_index', 'fs_delay', '-v7');
info = struct('impl', 'matlab', 'release', version, ...
    'commit', git_short(getenv_or(fullfile(here, '..', 'replay_matlab'), 'UWA_MATLAB_REPO')), ...
    'resamp_p', p, 'resamp_q', q);
save(fullfile(artifacts, 'matlab_replay.mat'), 'info', '-append');

%% Noise
noise = load(noise_path);
idx = cfg.noise.array_index(:).' + 1;
rng(cfg.noise.seed, 'twister');
w = noisegen([cfg.noise.n_samples, numel(idx)], cfg.noise.fs, idx, noise);
beta = noise.beta;

fprintf('matlab noisegen: w %dx%d, beta %dx%dx%d\n', size(w, 1), size(w, 2), size(beta));

Fs = noise.Fs; alpha = noise.alpha; fs_noise = cfg.noise.fs;
save(fullfile(artifacts, 'matlab_noise.mat'), 'w', 'beta', 'Fs', 'alpha', ...
    'fs_noise', 'idx', '-v7');
end


function v = getenv_or(default_value, name)
v = getenv(name);
if isempty(v); v = default_value; end
end


function sha = git_short(repo)
sha = 'unknown';
try
    [status, out] = system(sprintf('git -C "%s" rev-parse --short HEAD', repo));
    if status == 0; sha = strtrim(out); end
catch
end
end

% [EOF]
