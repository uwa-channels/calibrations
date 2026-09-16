function run_matlab()
% RUN_MATLAB Run the MATLAB implementation and write its results to artifacts/.
%
% Reads the same config.json and the same probe.mat as the Python runner, so
% neither language owns a second copy of the parameters or draws its own
% random numbers for the probe.
%
% For each case in config.json, produces
%    artifacts/matlab_<case>_replay.mat   y, y_same_start, y_static, y_resamp
%    artifacts/matlab_<case>_unpack.mat   u, u_fr
%    artifacts/matlab_<case>_noise.mat    w, beta
%
% `start` and `array_index` in config.json are the Python 0-based values, so
% MATLAB is given them plus one throughout.  A second replay is run at the same
% integer, which misaligns the two by a sample of fs_delay on purpose; the
% comparer checks that it looks clearly worse, which is how the harness
% establishes that it can resolve a difference that small.
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

probe = load(fullfile(artifacts, 'probe.mat'));
x = double(probe.input(:));

info = struct('impl', 'matlab', 'release', version, 'commit', ...
    git_short(getenv_or(fullfile(here, '..', 'replay_matlab'), 'UWA_MATLAB_REPO')));
save(fullfile(artifacts, 'matlab_env.mat'), 'info', '-v7');

want = strtrim(getenv('UWA_CALIBRATION_CASES'));
for c = 1:numel(cfg.data.cases)
    this = case_at(cfg.data.cases, c);
    if ~isempty(want) && ~ismember(this.name, strtrim(strsplit(want, ',')))
        continue
    end
    run_one(this, cfg, x, data, artifacts);
end
end


function c = case_at(cases, i)
% jsondecode returns a struct array when every case carries the same fields
% and a cell array when they do not, so accept either.
if iscell(cases)
    c = cases{i};
else
    c = cases(i);
end
end


function run_one(case_cfg, cfg, x, data, artifacts)
name = case_cfg.name;
channel = load(fullfile(data, case_cfg.channel_file));
noise = load(fullfile(data, case_cfg.noise_file));

fs = cfg.probe.fs;
fs_delay = channel.params.fs_delay;

%% Replay
array_index = case_cfg.replay.array_index(:).' + 1;   % config is 0-based
start = case_cfg.replay.start;

[p, q] = rat(fs_delay/fs);
y_resamp = resample(resample(x, p, q), q, p);

y = replay(x, fs, array_index, channel, start+1);
y_same_start = replay(x, fs, array_index, channel, start);

% The same replay with the phase trajectory removed, so that h_hat, the spline
% interpolation and the two resamplings are the only things acting.
static = channel;
for f = {'phi_hat', 'theta_hat', 'meta'}
    if isfield(static, f{1}); static = rmfield(static, f{1}); end
end
y_static = replay(x, fs, array_index, static, start+1);

fprintf('[%s] replay: y %dx%d, control %d, p/q = %d/%d\n', ...
    name, size(y, 1), size(y, 2), numel(y_resamp), p, q);
save(fullfile(artifacts, ['matlab_' name '_replay.mat']), 'y', ...
    'y_same_start', 'y_static', 'y_resamp', 'start', 'array_index', ...
    'fs_delay', '-v7');

%% Unpack, without and with f_resamp
uidx = case_cfg.unpack.array_index(:).' + 1;
tracked = channel;
if isfield(tracked, 'meta'); tracked = rmfield(tracked, 'meta'); end
u = unpack(case_cfg.unpack.fs_out, uidx, tracked, ...
    cfg.unpack.buffer_left, cfg.unpack.buffer_right);
tracked.f_resamp = case_cfg.unpack.f_resamp;
u_fr = unpack(case_cfg.unpack.fs_out, uidx, tracked, ...
    cfg.unpack.buffer_left, cfg.unpack.buffer_right);

fprintf('[%s] unpack: %dx%dx%d at %g Hz, f_resamp %g\n', name, size(u), ...
    case_cfg.unpack.fs_out, case_cfg.unpack.f_resamp);
save(fullfile(artifacts, ['matlab_' name '_unpack.mat']), 'u', 'u_fr', '-v7.3');

%% Noise
nidx = case_cfg.noise.array_index(:).' + 1;
rng(cfg.noise.seed, 'twister');
w = noisegen([cfg.noise.n_samples, numel(nidx)], cfg.noise.fs, nidx, noise);
beta = noise.beta;
fprintf('[%s] noisegen: w %dx%d, beta %dx%dx%d\n', name, size(w), size(beta));
save(fullfile(artifacts, ['matlab_' name '_noise.mat']), 'w', 'beta', '-v7');
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
