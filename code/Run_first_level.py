# ======================================================================
import os
import numpy as np
from spm import (
    Array, Cell, Struct, Runtime,
    spm_select, spm_dcm_specify, spm_dcm_load, spm_dcm_fit, spm_dcm_fmri_check
)

num = np.asarray

# Save the field of a (scalar) struct in a .mat file
_save = Runtime.call("eval", "@(f,x) save(f, '-struct', 'x')")
def save(*a): return _save(*a, nargout=0)

# Load a .mat into a (scalar) struct
def load(*a): return Runtime.call("load", *a)

# Change directory in MATLAB
def cd(*a): return Runtime.call("cd", *a)


# ======================================================================

this_dir = start_dir = os.path.abspath(os.path.dirname(__file__))

## Settings

# MRI scanner settings
TR = 3.6   # Repetition time (secs)
TE = 0.05  # Echo time (secs)

# Experiment settings
nsubjects   = 60
nregions    = 4
nconditions = 3

# Index of each condition in the DCM
# NOTE: Python uses 0-indexing
TASK, PICTURES, WORDS = 0, 1, 2

# Index of each region in the DCM
# NOTE: Python uses 0-indexing
lvF, ldF, rvF, rdF = 0, 1, 2, 3

## Specify DCMs (one per subject)

# A-matrix (on / off)
a = np.ones([nregions, nregions])
a[lvF, rdF] = 0
a[rdF, lvF] = 0
a[ldF, rvF] = 0
a[rvF, ldF] = 0

# B-matrix
b = np.zeros([nregions, nregions, 3])
b[:, :, TASK]     = np.zeros(nregions)  # Task
b[:, :, PICTURES] = np.eye(nregions)    # Pictures
b[:, :, WORDS]    = np.eye(nregions)    # Words

# C-matrix
c = np.zeros([nregions, nconditions])
c[:, TASK] = 1

# D-matrix (disabled)
d = np.zeros([nregions, nregions, 0])


for subject in range(nsubjects):

    name = f'sub-{subject+1:02d}'

    # Load SPM
    glm_dir = os.path.join(this_dir, '..', 'GLM', name)
    SPM     = load(os.path.join(glm_dir, 'SPM.mat'))
    SPM     = SPM.SPM

    # Load ROIs
    xY = Struct()
    f = [os.path.join(glm_dir, 'VOI_lvF_1.mat'),
         os.path.join(glm_dir, 'VOI_ldF_1.mat'),
         os.path.join(glm_dir, 'VOI_rvF_1.mat'),
         os.path.join(glm_dir, 'VOI_rdF_1.mat')]
    for r in range(len(f)):
        XY = load(f[r])
        xY[r] = XY.xY

        # Fix -- otherwise spm_dcm_specify crashes
        Ic = int(XY.xY.Ic.item())
        SPM.xCon[Ic].c = num([])

    # Move to output directory
    cd(glm_dir)

    # Select whether to include each condition from the design matrix
    # (Task, Pictures, Words)
    include = np.ones([3, 1])

    # Specify. Corresponds to the series of questions in the GUI.
    s = Struct()
    s.name       = 'full'
    s.u          = include              # Conditions
    s.delays     = num([TR]*nregions)   # Slice timing for each region
    s.TE         = TE
    s.nonlinear  = False
    s.two_state  = False
    s.stochastic = False
    s.centre     = True
    s.induced    = 0
    s.a          = a
    s.b          = b
    s.c          = c
    s.d          = d
    DCM = spm_dcm_specify(SPM, xY, s)

    # Return to script directory
    cd(start_dir)


## Collate into a GCM file and estimate

# Find all DCM files
dcms = spm_select('FPListRec', os.path.join(this_dir, '..', 'GLM'), 'DCM_full.mat')
dcms = dcms[:1]

# Prepare output directory
out_dir = os.path.join(this_dir, '..', 'analyses')
os.makedirs(out_dir, exist_ok=True)

# Check if it exists
if os.path.exists(os.path.join(out_dir, 'GCM_full.mat')):
    f = Runtime.call(
        'questdlg',
        'Overwrite existing GCM?', 'Overwrite?', 'Yes', 'No',
        {'Default': 'No', 'Interpreter': 'none'}
    )
    tf = (f == "Yes")
else:
    tf = True

# Collate & estimate
if tf:
    # Character array -> cell array
    # NOTE: In spm-python, `spm_select` already returns a cell array
    GCM = dcms

    # Filenames -> DCM structures
    GCM = spm_dcm_load(GCM)

    # Estimate DCMs (this won't effect original DCM files)
    GCM = spm_dcm_fit(GCM)

    # Save estimated GCM
    save(os.path.join(out_dir, 'GCM_full.mat'), {'GCM': GCM})


## Specify 28 alternative models structures
#  These will be templates for the group analysis

# Define B-matrix for each family (factor: task)
# -------------------------------------------------------------------------
# Both
b_task_fam = Cell()
b_task_fam[0] = np.zeros([4, 4, 2])
b_task_fam[0][:, :, 0] = 1   # Objects
b_task_fam[0][:, :, 1] = 1   # Words

# Words
b_task_fam[1] = np.zeros([4, 4, 2])
b_task_fam[1][:, :, 0] = 0   # Objects
b_task_fam[1][:, :, 1] = 1   # Words

# Objects
b_task_fam[2] = np.zeros([4, 4, 2])
b_task_fam[2][:, :, 0] = 1   # Objects
b_task_fam[2][:, :, 1] = 0   # Words

task_fam_names = ['Both', 'Words', 'Objects']

# Define B-matrix for each family (factor: dorsal-ventral)
# -------------------------------------------------------------------------
# Both
b_dv_fam = Cell()
b_dv_fam[0] = np.eye(4)

# Dorsal
b_dv_fam[1] = num([[0, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 1]])
# Ventral
b_dv_fam[2] = num([[1, 0, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 0]])

b_dv_fam_names = ['Both', 'Dorsal', 'Ventral']

# Define B-matrix for each family (factor: left-right)
# -------------------------------------------------------------------------
# Both
b_lr_fam = Cell()
b_lr_fam[0] = np.eye(4)

# Left
b_lr_fam[1] = num([[1, 0, 0, 0],
                   [0, 1, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 0, 0]])

# Right
b_lr_fam[2] = num([[0, 0, 0, 0],
                   [0, 0, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, 1]])

b_lr_fam_names = ['Both', 'Left', 'Right']

# Make a DCM for each mixture of these factors
# -------------------------------------------------------------------------

# Load and unpack an example DCM
GCM_full = load(os.path.join(out_dir, 'GCM_full.mat'))
GCM_full = spm_dcm_load(GCM_full.GCM)
DCM_template = GCM_full[0]
a = DCM_template.a
c = DCM_template.c
d = DCM_template.d
options = DCM_template.options

# Output cell array for new models
GCM_templates = Cell()

m = 0
task_family = Array()
b_dv_family = Array()
b_lr_family = Array()
for t in range(len(b_task_fam)):
    for dv in range(len(b_dv_fam)):
        for lr in range(len(b_lr_fam)):

            # Prepare B-matrix
            b = np.zeros([4, 4, 3])
            b[:, :, 1:] = (b_dv_fam[dv][..., None] + b_lr_fam[lr][..., None] + b_task_fam[t]) > 0

            # Prepare model name
            name = (
                'Task: %s, Dorsoventral: %s, Hemi: %s' %
                (task_fam_names[t], b_dv_fam_names[dv], b_lr_fam_names[lr])
            )

            # Build minimal DCM
            DCM = Struct()
            DCM.a       = a
            DCM.b       = b
            DCM.c       = c
            DCM.d       = d
            DCM.options = options
            DCM.name    = name
            GCM_templates[m] = DCM

            # Record the assignment of this model to each family
            task_family[m] = t
            b_dv_family[m] = dv
            b_lr_family[m] = lr
            m += 1


# Add a null model with no modulation
# -------------------------------------------------------------------------
b = np.zeros(4)
c = num([[1, 0, 0],
         [1, 0, 0],
         [1, 0, 0],
         [1, 0, 0]])
name = 'Task: None'

DCM = Struct()
DCM.b = np.zeros([4, 4, 3])
DCM.b[:, :, 1] = b
DCM.b[:, :, 2] = b
DCM.c          = c
DCM.name       = name

GCM_templates[m] = DCM

# Record the assignment of this model to each family
b_dv_family[m] = len(b_dv_fam)
b_lr_family[m] = len(b_lr_fam)
task_family[m] = len(b_task_fam)

m += 1

# Save
GCM = GCM_templates
save(os.path.join(out_dir, 'GCM_templates.mat'),
     {
        'GCM': GCM,
        'task_family': task_family,
        'b_dv_family': b_dv_family,
        'b_lr_family': b_lr_family,
     }
)

## Run diagnostics
GCM = load(os.path.join(out_dir, 'GCM_full.mat')).GCM
spm_dcm_fmri_check(GCM)
