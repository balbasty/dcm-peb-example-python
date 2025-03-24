## Reproduces the figures from the DCM/PEB tutorial papers
#
# Please run the analyses first, by uncommenting the following 4 lines:
#
# Run_first_level;
# Run_second_level;
# close all;
# clear all;
#
# Then run the code below.

import numpy as np
from spm import Runtime, spm_vec, spm_find_pC, spm_plot_ci, spm_dcm_fmri_check
from pz_utils import pz_dcm_get_parameter_names, pz_rescale


class Builtin:
    # Wrapper for calls to builtin functions

    class WrappedBuiltin:

        def __init__(self, name):
            self.name = name

        def __call__(self, *args, **kwargs):
            return Runtime.call(self.name, *args, **kwargs)

    def __getattr__(self, name):
        return self.WrappedBuiltin(name)

    def import_in(self, env, *args, **kwargs):
        for name in args:
            env[name] = getattr(self, name)
        for as_name, name in kwargs.items():
            env[as_name] = getattr(self, name)


plt = builtin = Builtin()

# Save the field of a (scalar) struct in a .mat file
save = builtin.eval("@(f,x) save(f,'-struct',x)")

# Load a .mat into a (scalar) struct
load = builtin.load

# Short alias to specify that something is a numeric array (and not a cell)
num = np.asarray

# --------------

## Load GCM
GCM = load('../analyses/GCM_full_pre_estimated.mat').GCM
X = load('../design_matrix.mat', 'X').X
## Show inputs matrix U and timeseries matrix Y. (Part 1, Figure 2)

DCM = GCM[0]

plt.figure(
    'Name', 'Input matrix U',
    'NumberTitle', 'off'
)
x = num(range(1, len(DCM.U.u)+1)) * DCM.U.dt
plt.imagesc(DCM.U.u > 0)
plt.colormap('gray')
plt.set(plt.gca(), 'YTickLabel', x[::500].round(), 'YTick', num(range(500))+1)
plt.set(plt.gca(), 'XTickLabel', ['Task','Pictures','Words'], 'XTick', num(range(3))+1)
plt.set(plt.gca(), 'FontSize', 12)
plt.ylabel('Time (secs)')

plt.figure('Name', 'Timeseries Y', 'NumberTitle', 'off')
x   = num(range(1, DCM.v+1)) * DCM.Y.dt
y   = DCM.Y.y
for i in range(1, 5):
    plt.subplot(1, 4, i)
    plt.plot(y[:, i], x, 'Color', 'k')
    plt.set(plt.gca(), 'Ydir', 'reverse')
    plt.set(plt.gca(), 'XLim', num([-2, 2]), 'XTick', num([]), 'FontSize', 12)
    plt.ylim(num([0, 750]))
    if i > 1:
        plt.set(plt.gca(),'YTick',num([]))

## Switched on vs switched off connection (Part 1, Figure 5)
plt.figure()
x = np.linspace(-4, 4, 800)

plt.subplot(1, 2, 1)
y = builtin.normpdf(x, 0, 1)
plt.area(x, y)
plt.axis('square')
plt.set(plt.gca(), 'YTick', num([]), 'FontSize', 12)
plt.title('Prior - Switched on')
plt.xlabel('Parameter value')
plt.ylabel('Probability')

plt.subplot(1, 2, 2)
y = builtin.normpdf(x, 0,0.0001)
plt.area(x, y)
plt.axis('square')
plt.title('Prior - Switched off')
plt.xlabel('Parameter value')
plt.ylabel('Probability')
plt.set(plt.gca(), 'YTick', num([]), 'FontSize', 12)
## Parameters and timeseries from example subject (Part 1, Figure 6)

# Subject
s = 37

# Unpack DCM
DCM = GCM[s]
Ep = spm_vec(DCM.Ep)
Vp = spm_vec(DCM.Vp)
pnames = pz_dcm_get_parameter_names(DCM)

# Re-order as listed in paper
idx_a = num(range( 1, 17)).reshape(DCM.a.shape)
idx_b = num(range(17, 65)).reshape(DCM.b.shape)
idx_c = num(range(65, 76)).reshape(DCM.c.shape)
idx   = np.concat([
    np.diag(idx_a),                             # A intrinsic
    spm_vec(idx_a - np.diag(np.diag(idx_a))),   # A extrinsic
    np.diag(idx_b[:, :, 1]),                    # B pictures intrinsic
    np.diag(idx_b[:, :, 2]),                    # B words intrinsic
    idx_c[:,0]                                  # C Task driving
])
idx = idx[idx != 0]

# Further limit to free parameters
idx_on = spm_find_pC(DCM) - 1  # NOTE: 1-indexing to 0-indexing
idx_both = num([])
for i in range(len(idx)):
    if idx[i] in idx_on:
        idx_both.append(idx(i))
idx = idx_both

# Filter
Ep      = Ep[idx]
Vp      = Vp[idx]
pnames  = pnames[idx]

# Plot
plt.figure(
    'Name', 'Parameters and timeseries of example subject',
    'NumberTitle', 'off'
)
plt.subplot(2, 1, 1)
Cp = np.diag(Vp)
spm_plot_ci(Ep, Cp)
plt.hold('on')
x = np.repeat(range(1, 5), len(Ep)/4)
plt.set(plt.gca(), 'XTick', num(range(len(Ep))) + 1, 'XTickLabel', x)

# Decorate
plt.xlabel('DCM parameter')
plt.ylabel('Posterior')

# Plot example timeseries
plt.subplot(2, 1, 2)
x = (1 + np.arange(DCM.v)) * DCM.Y.dt
plt.plot(x, DCM.y, 'LineWidth', 2)
plt.hold('on')
plt.plot(x,DCM.y + DCM.R[:, :4], ':')
plt.legend(['lvF', 'ldF', 'rvF', 'rdF'])
plt.hold('on')
x = (1 + np.arange(DCM.U.u)) * DCM.U.dt
plt.plot(x, builtin.full(DCM.U.u[:, 0]), 'Color', num([0.2, 0.2, 0.2]))
## Report explained variance across subjects
GCM_diagnostics = spm_dcm_fmri_check(GCM)
exp_var = num([x.diagnostics(1) for x in GCM_diagnostics])
print('Mean explained variance: %2.2f std: %2.2f' % (exp_var.mean(), exp_var.std()))
## Covariance components (Part 1, Figure A.1)
plt.figure('Name', 'DCM Covariance components', 'NumberTitle', 'off')
for i in range(4):
    plt.subplot(1, 4, i+1);
    plt.imagesc(GCM[0].Y.Q[i])
    plt.axis('tight')
    plt.axis('square')
    plt.axis('off')
    plt.colormap('gray')
## Plot cartoon of T2* relaxation (Part 1, Figure A.3)
B0 = 3
gamma = 42.56        # Gyromagnetic ratio (Mhz/Tesla)
larmor = gamma * B0  # MHz
T2 = 0.2             # Transverse time constant, arbitrarily chosen for display (secs)
flipangle = np.pi/2  # Radians

# Envelope
t = np.linspace(0, 1, 1000)  # secs
envelope = np.sin(flipangle) * np.exp(-t / T2)

plt.figure()
plt.plot(t, envelope, 'Color', 'k')
plt.hold('on')

# Sine wave
y = np.sin(flipangle) * np.sin(larmor * t) * np.exp(-t/T2)
plt.plot(t, y, 'Color', 'k', 'LineWidth', 2)
plt.xlabel('Time')
plt.ylabel('Measured signal')

plt.line(num([T2, T2]), num([envelope(t==T2), 1]), 'Color', 'k')
plt.set(plt.gca(), 'XTick', [], 'YTick', [])
## PEB design matrices in colour. (Part 2 Figure 3)

colours = num([
    [0, 0, 0],
    [100, 100, 100],

    [98, 62.4, 64.7],   # Reds: 2->6
    [92.9, 43.1, 46.7],
    [80, 26.7, 30.6],
    [68.2, 14.5, 18.4],
    [53.7, 4.3, 7.8],

    [100, 87.8, 63.1],  # Yellows: 7->11
    [95.3, 78.4, 44.3],
    [82, 63.5, 27.1],
    [70.2, 51.8, 14.9],
    [55.3, 38.4, 4.3],

    [51.4, 58.8, 76.5], # Blues: 12:
    [33.3, 42.4, 63.9],
    [22, 32.2, 54.9],
    [13.7, 23.9, 47.1],
    [6.7, 15.7, 36.9],

    [60, 87.5, 55.3], # greens
    [42.4, 78.8, 36.9],
    [28.6, 67.8, 22.4],
    [18.4, 58, 12.2],
    [9.4, 45.9, 3.5],
]) / 100

reds    = num(range(2, 8))
yellows = num(range(8, 13))
blues   = num(range(13, 18))
greens  = num(range(18, 23))

PEB_B = load('../analyses/PEB_B.mat').PEB_B

# We z-score for display purposes - it's not zscored in the analysis
PEB.M.X(:,2:end) = zscore(PEB.M.X(:,2:end) );

XB = PEB.M.X;
XW = PEB.M.W;

XB(:,2) = pz_rescale(XB(:,2),reds(1),reds(end)-0.01);
XB(:,3) = pz_rescale(XB(:,3),yellows(1),yellows(end)-0.01);
XB(:,4) = pz_rescale(XB(:,4),blues(1),blues(end)-0.01);
XB(:,5) = pz_rescale(XB(:,5),greens(1),greens(end)-0.01);

X = kron(XB,XW);

XB(:,1) = 2;

figure('Name','PEB design matrices','NumberTitle','off');
subplot(1,3,1);
image(XB); colormap(gca,colours); axis square;
xlabel('Covariate'); ylabel('Subject');
set(gca,'FontSize',12);
title('Between-Subjects X_B','FontSize',16);

subplot(1,3,2);
imagesc(XW); colormap gray; axis square;
xlabel('DCM parameter'); ylabel('DCM parameter');
set(gca,'FontSize',12);
title('Within-Subjects X_W','FontSize',16);

subplot(1,3,3);
imagesc(X);  colormap(gca,colours); axis square;
xlabel('Group level covariate'); ylabel('Subject level DCM parameter');
set(gca,'FontSize',12);
title('Design matrix X','FontSize',16);

## PEB GLM parameters (Part 2, Figure 4)
load('../analyses/PEB_B.mat','PEB_B');

nconnections = length(PEB_B.Pnames);

# One bar plot for mean and LI
nx = 2;
Ep = PEB_B.Ep(:,1:nx);
Ep = Ep(:);
Vp = diag(PEB_B.Cp);
Vp = Vp(1:(nconnections * nx));

figure('Name','PEB GLM parameters','NumberTitle','off');

subplot(1,2,1);
spm_plot_ci(Ep,diag(Vp));
xlabel('GLM Parameter');
ylabel('Estimate');
set(gca,'FontSize',12);
hold on;
x=1:nconnections:(nconnections*nx);
for i = 2:length(x)
    line([x(i)-0.5 x(i)-0.5],[-2 4],'Color',[0.2 0.2 0.2]);
end
axis square;
title('Group level GLM parameters \theta^{(2)}');

subplot(1,2,2);
bar(diag(PEB_B.Ce),'FaceColor',[1 1 1]*.8);
axis square;
xlabel('Connectivity Parameter');
ylabel('Estimate');
set(gca,'FontSize',12);
title('Random effects variance diag(\Sigma^{(2)})');
## Explicit model comparison figure (Part 2, Figure 5)
load('../analyses/BMA_B_28models.mat','BMA');

figure('Name','Comparison of pre-defined models','NumberTitle','off');

# Plot model space
subplot(2,3,1);
imagesc(BMA.K)
axis square; colormap gray;
set(gca,'XTick',[],'YTick',[]);

Pp_common     = sum(BMA.P,2);
Pp_diff = sum(BMA.P,1);

subplot(2,3,2);
imagesc(BMA.P);
axis square;
xlabel('Model (differences)'); ylabel('Model (commonalities)');
title('Posterior probabilities','FontSize',16);
set(gca,'FontSize',12);
colormap gray;

subplot(2,3,3);
bar(Pp_common); ylim([0 1]); xlim([0 size(Pp_common,1)]);
axis square;
xlabel('Model'); ylabel('Probability');
title('Commonalities','FontSize',16);
set(gca,'FontSize',12);

subplot(2,3,4);
bar(Pp_diff);  ylim([0 1]); xlim([0 size(Pp_common,1)]);
axis square;
xlabel('Model');  ylabel('Probability');
title('Differences (LI)','FontSize',16);
set(gca,'FontSize',12);

# Show connections in winning model 4
BMA.Kname(BMA.K(4,:)==1)

# Show connections in winning model 15
BMA.Kname(BMA.K(15,:)==1)

# BMA (specific models) before and after thresholding
nconnections = length(BMA.Pnames);

# One bar plot for mean and LI
nx = 2;
Ep = BMA.Ep(1:(nconnections * nx));
Ep = Ep(:);
Vp = BMA.Cp;
Vp = Vp(1:(nconnections * nx));

subplot(2,3,5);
spm_plot_ci(Ep,diag(Vp));
xlabel('GLM Parameter');
ylabel('Estimate');
set(gca,'FontSize',12);
hold on;
x=1:nconnections:(nconnections*nx);
for i = 2:length(x)
    line([x(i)-0.5 x(i)-0.5],[-2 4],'Color',[0.2 0.2 0.2]);
end
axis square;
title('Bayesian Model Average');

# Threshold commonalities
Ep(1:nconnections) = Ep(1:nconnections) .* (BMA.Pw > 0.95)';
Vp(1:nconnections) = Vp(1:nconnections) .* (BMA.Pw > 0.95)';

# Threshold differences
Ep(nconnections+1:end) = Ep(nconnections+1:end) .* (BMA.Px > 0.95)';
Vp(nconnections+1:end) = Vp(nconnections+1:end) .* (BMA.Px > 0.95)';

subplot(2,3,6);
spm_plot_ci(Ep,diag(Vp));
xlabel('GLM Parameter');
ylabel('Estimate');
set(gca,'FontSize',12);
hold on;
x=1:nconnections:(nconnections*nx);
for i = 2:length(x)
    line([x(i)-0.5 x(i)-0.5],[-2 4],'Color',[0.2 0.2 0.2]);
end
axis square;
title('Thresholded');
## Family comparisons (Part 2, Figure 6)
load('../analyses/BMA_fam_task.mat');
load('../analyses/BMA_fam_b_dv.mat');
load('../analyses/BMA_fam_b_lr.mat');

rows = 1; cols = 3;

figure('Name','Family comparison','NumberTitle','off');

subplot(rows,cols,1);
imagesc(fam_task.family.post); axis square; colormap gray;
ylabel('Commonalities');xlabel('Differences');
title('Factor 1: Task','FontSize',16);
set(gca,'FontSize',12,'YTick',1:size(fam_task.family.post,1));

subplot(rows,cols,2);
imagesc(fam_b_dv.family.post); axis square; colormap gray;
ylabel('Commonalities');xlabel('Differences');
title('Factor 2: Dorsoventral','FontSize',16);
set(gca,'FontSize',12,'YTick',1:size(fam_b_dv.family.post,1));

subplot(rows,cols,3);
imagesc(fam_b_lr.family.post); axis square; colormap gray;
ylabel('Commonalities');xlabel('Differences');
title('Factor 3: Left/right','FontSize',16);
set(gca,'FontSize',12,'YTick',1:size(fam_b_lr.family.post,1));
## BMA after automatic search (Part 2, Figure 7)

load('../analyses/BMA_search_B.mat','BMA_B');

nconnections = length(BMA_B.Pnames);

# One bar plot for mean and LI
nx = 2;
p  = 1:(nconnections * nx);
Ep = BMA_B.Ep(p);
Pp = BMA_B.Pp(p);
Vp = diag(BMA_B.Cp);
Vp = Vp(p);

figure('Name','BMA - automatic search','NumberTitle','off');

subplot(1,2,1);
spm_plot_ci(Ep,diag(Vp));
xlabel('GLM Parameter');
ylabel('Estimate');
set(gca,'FontSize',12);
hold on;
x=1:nconnections:(nconnections*nx);
for i = 2:length(x)
    line([x(i)-0.5 x(i)-0.5],[-2 4],'Color',[0.2 0.2 0.2]);
end
axis square;
title('Bayesian Model Average');

# Threshold
Ep = Ep .* (Pp(:) > 0.95);
Vp = Vp .* (Pp(:) > 0.95);

subplot(1,2,2);
spm_plot_ci(Ep,diag(Vp));
xlabel('GLM Parameter');
ylabel('Estimate');
set(gca,'FontSize',12);
hold on;
x=1:nconnections:(nconnections*nx);
for i = 2:length(x)
    line([x(i)-0.5 x(i)-0.5],[-2 4],'Color',[0.2 0.2 0.2]);
end
axis square;
title('Thresholded');

# Show connections in BMA
BMA.Kname(BMA.K(4,:)==1)


## Leave-one-out cross validation (Part 2, Figure 8)
load('../analyses/LOO_rdF_words.mat','qE','qC','Q');
load('../design_matrix.mat','X');

# Subject order
k = 1:size(X,1);

figure('Name','Leave one out cross validation)','NumberTitle','off');

figure;
subplot(1,2,1), spm_plot_ci(qE(k),qC(k)), hold on
plot(X(k,2),'--','Color','K','LineWidth',2), hold off
xlabel('Subject','FontSize',12), ylabel('Predicted subject effect','FontSize',12)
title('Out of sample estimates','FontSize',16)
axis tight, axis square;
set(gca,'FontSize',12);

# count the number of subjects where preditcion was within 99# CI
#--------------------------------------------------------------------------
ci = 0.99;
ci = 1 - (1-ci)/2;
ci = spm_invNcdf(ci);
c  = ci*sqrt(qC);
lower  = qE - c;
higher = qE + c;

sum(X(:,2)' >= lower & X(:,2)' <= higher)

# classical inference on classification accuracy
#--------------------------------------------------------------------------
[T,df] = spm_ancova(X(:,1:2),[],qE(:),[0;1]);
r      = corrcoef(qE(:),X(:,2));
r      = full(r(1,2));

if isnan(T)
    p = NaN;
else
    p = 1 - spm_Tcdf(T,df(2));
end
str = sprintf('corr(df:#-2.0f) = #-0.2f: p = #-0.5f',df(2),r,p);

subplot(1,2,2)
plot(X(:,2),qE,'.','Markersize',8,'Color','k')
xlabel('Group effect','FontSize',12), ylabel('Estimate','FontSize',12)
title(str,'FontSize',16)
set(gca,'FontSize',12);
axis square;lsline;

## PEB precision components (Part 2, Figure 9)
load('../analyses/PEB_B.mat');
PEB = PEB_B;

n = 4;

nq = length(PEB.M.Q);
ns = size(GCM,1);
np = nq; # number of dcm parameters

figure;

subplot(1,n+2,1);
imagesc( kron(eye(ns), np) );
colormap gray;
axis square;
set(gca,'XTick',[],'YTick',[]);

subplot(1,n+2,2);
Q0=eye(np);
imagesc(Q0);
colormap gray;
axis square;
set(gca,'XTick',[],'YTick',[]);

i = 3;
for q = [1:n-1 nq]
    Q = PEB.M.Q{q};
    subplot(1,n+2,i);
    imagesc(Q);
    colormap gray;
    axis square;
    set(gca,'XTick',[],'YTick',[]);
    i = i + 1;
end
