%% THORPE LENGTH
% Calculating the Lt (Thorpe lengthscale) for a given DZ.
% Note we search for the displacement using the original profile (not binned).
% CT = original data set
% CT_dz = CT
% CT_dz_sort = sorted

tmp = importdata("data\L0\t_88-22-200.txt").data;

% crn = tmp(:,1);
% day = tmp(:,2);
% p = tmp(:,3);
T = tmp(:,4);

Lt      = NaN(zl,1);
for n=1:zl
    [~,idx] = min(abs(CT_dz_sort(n) - CT));
    Lt(n)   = abs(zb(n) - p(idx));
    clear idx
    % should be p(n) - p(idx)
end

% SORTED PROFILES
CT_dz_sort  = sort(CT_dz,'descend');

% Now calculate epsilon and diffusion
eps_thorpe      = 0.640 .* Lt.^2 .* sqrt(N2_dzs).^3;
D_thorpe        = 0.128 .* Lt.^2 .* sqrt(N2_dzs);