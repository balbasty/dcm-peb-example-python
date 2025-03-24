# ======================================================================
import os
import shutil
from spm import (
    Struct,
    Runtime,
    spm_dcm_peb,
    spm_dcm_peb_bmc,
    spm_dcm_peb_bmc_fam,
    spm_dcm_loo,
    spm_dcm_peb_review,
)


# Save the field of a (scalar) struct in a .mat file
def save(f, x): return Runtime.call("save", f, "-struct", x)


# Load a .mat into a (scalar) struct
def load(*a): return Runtime.call("load", *a)


# ======================================================================

# ----------------------
# Load PEB prerequisites
# ----------------------

# Load design matrix
dm = load('../design_matrix.mat')
X = dm.X
X_labels = dm.labels

# Import downloaded GCM file if needed
if not os.path.exists('../analyses/GCM_full.mat'):
    shutil.copyfile(
        '../analyses/GCM_full_pre_estimated.mat',
        '../analyses/GCM_full.mat'
    )

# Load GCM
GCM = load('../analyses/GCM_full.mat').GCM

# PEB settings
M = Struct()
M.Q = 'all'
M.X = X
M.Xnames = X_labels
M.maxit = 256

# ------------------------------
# Build PEB (using B parameters)
# ------------------------------

[PEB_B, RCM_B] = spm_dcm_peb(GCM, M, ['B'], nargout=2)
save('../analyses/PEB_B.mat', {'PEB_B': PEB_B, 'RCM_B': RCM_B})

# ----------------
# Automatic search
# ----------------

BMA_B = spm_dcm_peb_bmc(PEB_B)
save('../analyses/BMA_search_B.mat', {'BMA_B': BMA_B})

# -----------------------------
# Hypothesis-based analysis (B)
# -----------------------------

# Load estimated PEB
PEB_B = load('../analyses/PEB_B.mat').PEB_B

# Load template models
templates = load('../analyses/GCM_templates.mat')

# Run model comparison
[BMA, BMR] = spm_dcm_peb_bmc(PEB_B, templates.GCM, nargout=2)

# Show connections in winning model 4
BMA.Kname(BMA.K[3, :] == 1)

# Show connections in winning model 15
BMA.Kname(BMA.K[14, :] == 1)

save('../analyses/BMA_B_28models.mat', {'BMA': BMA, 'BMR': BMR})

# ---------------
# Family analysis
# ---------------

# Load the result from the comparison of 28 reduced models
_ = load('../analyses/BMA_B_28models.mat')

# Compare families
[BMA_fam_task, fam_task] = spm_dcm_peb_bmc_fam(BMA, BMR, templates.task_family, 'ALL', nargout=2)

[BMA_fam_b_dv, fam_b_dv] = spm_dcm_peb_bmc_fam(BMA, BMR, templates.b_dv_family, 'NONE', nargout=2)

[BMA_fam_b_lr, fam_b_lr] = spm_dcm_peb_bmc_fam(BMA, BMR, templates.b_lr_family, 'NONE', nargout=2)

save('../analyses/BMA_fam_task.mat', {'BMA_fam_task': BMA_fam_task, 'fam_task': fam_task})
save('../analyses/BMA_fam_b_dv.mat', {'BMA_fam_b_dv': BMA_fam_b_dv, 'fam_b_dv': fam_b_dv})
save('../analyses/BMA_fam_b_lr.mat', {'BMA_fam_b_lr': BMA_fam_b_lr, 'fam_b_lr': fam_b_lr})

# LOO
[qE, qC, Q] = spm_dcm_loo(GCM, M, ['B(4,4,3)'], nargout=3)
save('../analyses/LOO_rdF_words.mat', {'qE': qE, 'qC': qC, 'Q': Q})

# Correlate rdF
B = [gcm.Ep.B(3, 3, 2) for gcm in GCM]
LI = X[:, 1]
Runtime.call('figure')
Runtime.call('scatter', LI, B)
Runtime.call('lsline')
[R, P] = Runtime.call('corrcoeff', LI, B, nargout=2)

# Build PEB (A)
[PEB_A, RCM_A] = spm_dcm_peb(GCM[:, 0], M, ['A'], nargout=2)
save('../analyses/PEB_A.mat', {'PEB_A': PEB_A, 'RCM_A': RCM_A})

# Search-based analysis (A)
PEB_A = load('../analyses/PEB_A.mat').PEB_A
BMA_A = spm_dcm_peb_bmc(PEB_A)
save('../analyses/BMA_search_A.mat', {'BMA_A': BMA_A})
spm_dcm_peb_review(BMA_A, GCM)
