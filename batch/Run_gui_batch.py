# The script runs the step-by-step GUI example from start to end using the
# batch complete_gui_batch.mat.
# -------------------------------------------------------------------------
import os.path as op
from spm import spm_select, spm_jobman, Runtime

# Load a .mat into a (scalar) struct
load = lambda *a: Runtime.call("load", *a)


this_dir = op.abspath(op.dirname(__file__))
glm_dir = op.join(this_dir, '..', 'GLM')


# Where to store analysis results
dir_out = op.join(this_dir, '..', 'analyses')

# Where to find manually created template DCMs
dcm_full = op.join(glm_dir, 'sub-01', 'DCM_full.mat')
dcm_alt = op.join(glm_dir, 'sub-01', 'DCM_no_ldF_modulation.mat')

# Find SPM.mat files for all subjects (timing)
spms = spm_select('FPListRec', glm_dir, 'SPM.mat')
assert len(spms) == 60

# Find timeseries for all subjects
lvF = spm_select('FPListRec', glm_dir, 'VOI_lvF')
ldF = spm_select('FPListRec', glm_dir, 'VOI_ldF')
rvF = spm_select('FPListRec', glm_dir, 'VOI_rvF')
rdF = spm_select('FPListRec', glm_dir, 'VOI_rdF')
assert len(lvF) == 60
assert len(ldF) == 60
assert len(rvF) == 60
assert len(rdF) == 60

# Prepare the batch
matlabbatch = load(op.join(this_dir, 'complete_gui_batch.mat')).matlabbatch

matlabbatch[0].spm.dcm.spec.fmri.group.output.dir = [dir_out]
matlabbatch[0].spm.dcm.spec.fmri.group.template.fulldcm = [dcm_full]
matlabbatch[0].spm.dcm.spec.fmri.group.template.altdcm = [dcm_alt]
matlabbatch[0].spm.dcm.spec.fmri.group.data.spmmats = spms
matlabbatch[0].spm.dcm.spec.fmri.group.data.region = [lvF, ldF, rvF, rdF]

# Run
spm_jobman('run', matlabbatch)
