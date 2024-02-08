
% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Wed Feb  7 19:10:54 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.68554E-02 0.00140  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73145E-01 3.9E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.33637E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98725E-01 7.1E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98208E-01 9.9E-06 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.33293E-01 0.00014  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.83747E+00 0.00053  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.70509E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.70509E+01 0.00051  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.79200E+00 0.00073  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.29936E-01 0.00155  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000421 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00084E+03 0.00124 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00084E+03 0.00124 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.00418E+01 ;
RUNNING_TIME              (idx, 1)        =  2.15080E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.20667E-02  2.20667E-02 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  1.96296E+01  1.96296E+01  0.00000E+00 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.13010E+01  0.00000E+00 ];
CPU_USAGE                 (idx, 1)        = 0.93183 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.97537E-01 0.00183 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  8.52608E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  2.20000E-04 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 23 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  8.20680E+05 ;
TOT_DECAY_HEAT            (idx, 1)        =  6.31024E-07 ;
TOT_SF_RATE               (idx, 1)        =  4.05567E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  8.20680E+05 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  6.31024E-07 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  0.00000E+00 ;
INHALATION_TOXICITY       (idx, 1)        =  7.58853E+00 ;
INGESTION_TOXICITY        (idx, 1)        =  3.98676E-02 ;
ACTINIDE_INH_TOX          (idx, 1)        =  7.58853E+00 ;
ACTINIDE_ING_TOX          (idx, 1)        =  3.98676E-02 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  0.00000E+00 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  0.00000E+00 ;
SR90_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
TE132_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
I131_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
I132_ACTIVITY             (idx, 1)        =  0.00000E+00 ;
CS134_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
CS137_ACTIVITY            (idx, 1)        =  0.00000E+00 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  1.19575E+05 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.11108E-02 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  8.20602E+05 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  6.39276E+05 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.67433E+09 0.00094  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 0 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
BURN_DAYS                 (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
FIMA                      (idx, [1:  3])  = [  0.00000E+00  0.00000E+00  1.58926E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.84006E-01 0.00229 ];
U235_FISS                 (idx, [1:   4]) = [  6.16099E+12 0.00100  9.46812E-01 0.00030 ];
U238_FISS                 (idx, [1:   4]) = [  3.45456E+11 0.00574  5.30510E-02 0.00541 ];
U235_CAPT                 (idx, [1:   4]) = [  1.46064E+12 0.00256  3.00719E-01 0.00222 ];
U238_CAPT                 (idx, [1:   4]) = [  2.89150E+12 0.00215  5.95088E-01 0.00122 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000421 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51301E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000421 1.00151E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 427566 4.28040E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 572855 5.73473E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000421 1.00151E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.04774E-09 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.59722E+13 2.3E-05  1.59722E+13 2.3E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49727E+12 1.9E-06  6.49727E+12 1.9E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  4.86133E+12 0.00096  4.40619E+12 0.00104  4.55149E+11 0.00117 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.13586E+13 0.00041  1.09035E+13 0.00042  4.55149E+11 0.00117 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.13487E+13 0.00094  1.13487E+13 0.00094  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.79115E+14 0.00082  1.59147E+14 0.00090  3.19968E+14 0.00085 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.13586E+13 0.00041 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.06996E+14 0.00067 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27837E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.90339E+00 0.00065 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.52528E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.77012E-01 0.00082 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.34806E+00 0.00080 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.40986E+00 0.00083 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.40986E+00 0.00083 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45829E+00 2.4E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02468E+02 1.9E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.40959E+00 0.00087  1.40010E+00 0.00084  9.76443E-03 0.01536 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.40844E+00 0.00041 ];
COL_KEFF                  (idx, [1:   2]) = [  1.40802E+00 0.00093 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.40844E+00 0.00041 ];
ABS_KINF                  (idx, [1:   2]) = [  1.40844E+00 0.00041 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.70514E+01 0.00039 ];
IMP_ALF                   (idx, [1:   2]) = [  1.70559E+01 0.00017 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  7.95163E-07 0.00665 ];
IMP_EALF                  (idx, [1:   2]) = [  7.84633E-07 0.00294 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.90184E-01 0.00553 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.88930E-01 0.00243 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.04056E-03 0.01156  1.59052E-04 0.07124  6.91034E-04 0.03324  4.82248E-04 0.03914  9.97522E-04 0.02621  1.59704E-03 0.02059  5.11615E-04 0.03509  4.29553E-04 0.04122  1.72493E-04 0.06168 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.55790E-01 0.01841  4.21374E-03 0.06265  2.37650E-02 0.01954  3.03624E-02 0.02833  1.25592E-01 0.01090  2.90712E-01 0.00348  5.23859E-01 0.02336  1.15089E+00 0.02903  1.40762E+00 0.05529 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  6.99831E-03 0.01642  2.04040E-04 0.10772  9.40448E-04 0.04413  7.03475E-04 0.05497  1.36377E-03 0.03634  2.22262E-03 0.02795  7.08725E-04 0.05045  6.05005E-04 0.05502  2.50219E-04 0.09532 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.56703E-01 0.02656  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.0E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.22675E-05 0.00209  1.22599E-05 0.00209  1.34095E-05 0.02201 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.72854E-05 0.00188  1.72748E-05 0.00189  1.88809E-05 0.02193 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.90859E-03 0.01605  1.85781E-04 0.10378  9.77040E-04 0.04456  6.65361E-04 0.05236  1.37381E-03 0.03427  2.20652E-03 0.02811  7.10291E-04 0.04968  5.44543E-04 0.05896  2.45237E-04 0.08528 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.47068E-01 0.02728  1.24667E-02 0.0E+00  2.82917E-02 2.6E-09  4.25244E-02 7.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.6E-09  3.55460E+00 4.8E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.22254E-05 0.00447  1.22253E-05 0.00449  1.03354E-05 0.05191 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.72245E-05 0.00434  1.72244E-05 0.00435  1.45721E-05 0.05189 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  6.32376E-03 0.04479  8.30805E-05 0.41566  1.01567E-03 0.11625  6.77498E-04 0.13419  1.25827E-03 0.10396  1.92899E-03 0.08811  6.45347E-04 0.14659  3.81482E-04 0.17184  3.33418E-04 0.27074 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.42446E-01 0.07207  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.9E-09  2.92467E-01 5.8E-09  6.66488E-01 5.9E-09  1.63478E+00 0.0E+00  3.55460E+00 5.4E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  6.36509E-03 0.04324  7.69605E-05 0.43769  9.96020E-04 0.11384  7.06636E-04 0.13696  1.26414E-03 0.09826  1.92414E-03 0.08648  6.54695E-04 0.13838  4.16727E-04 0.16453  3.25775E-04 0.24678 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.48591E-01 0.07113  1.24667E-02 4.0E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.9E-09  2.92467E-01 5.8E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 3.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -5.23245E+02 0.04427 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.22943E-05 0.00120 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.73237E-05 0.00087 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.83306E-03 0.00785 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.56217E+02 0.00795 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.20395E-07 0.00115 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92419E-06 0.00094  2.92455E-06 0.00094  2.87061E-06 0.01222 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.99045E-05 0.00126  1.99035E-05 0.00127  2.01773E-05 0.01623 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.77689E-01 0.00082  5.76278E-01 0.00084  9.26704E-01 0.02147 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.26616E+01 0.02573 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.70164E+01 0.00051  2.97377E+01 0.00065 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.31051E+04 0.00783  5.37913E+04 0.00294  1.14198E+05 0.00259  1.25664E+05 0.00144  1.18496E+05 0.00194  1.33171E+05 0.00161  9.01691E+04 0.00137  8.10232E+04 0.00145  6.17020E+04 0.00155  5.01567E+04 0.00131  4.31400E+04 0.00245  3.91984E+04 0.00163  3.58646E+04 0.00191  3.41315E+04 0.00213  3.31781E+04 0.00192  2.85833E+04 0.00173  2.82450E+04 0.00258  2.78159E+04 0.00235  2.72157E+04 0.00124  5.26545E+04 0.00198  5.01378E+04 0.00158  3.57022E+04 0.00169  2.27908E+04 0.00251  2.58184E+04 0.00226  2.39805E+04 0.00185  2.17960E+04 0.00218  3.48008E+04 0.00160  8.04351E+03 0.00364  1.00426E+04 0.00345  9.27797E+03 0.00412  5.30996E+03 0.00424  9.30855E+03 0.00481  6.27651E+03 0.00478  5.23278E+03 0.00545  9.79145E+02 0.00841  9.87951E+02 0.00929  1.00911E+03 0.00881  1.03853E+03 0.00947  1.03260E+03 0.00939  1.01174E+03 0.00978  1.05311E+03 0.00787  9.82603E+02 0.00985  1.86507E+03 0.00739  3.02898E+03 0.00500  3.83773E+03 0.00429  1.00101E+04 0.00366  1.04214E+04 0.00298  1.09606E+04 0.00325  6.93550E+03 0.00404  4.84933E+03 0.00404  3.63242E+03 0.00461  4.11481E+03 0.00452  7.32704E+03 0.00292  9.12365E+03 0.00300  1.64057E+04 0.00337  2.23693E+04 0.00211  3.00384E+04 0.00144  1.78941E+04 0.00277  1.21900E+04 0.00286  8.52801E+03 0.00319  7.43068E+03 0.00255  7.13372E+03 0.00352  5.84503E+03 0.00398  3.85795E+03 0.00474  3.54973E+03 0.00273  3.11377E+03 0.00424  2.57334E+03 0.00416  1.97946E+03 0.00444  1.29658E+03 0.00528  4.49685E+02 0.00895 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.40798E+00 0.00098 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.21845E+14 0.00080  5.73154E+13 0.00081 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.46171E-01 0.00018  1.33727E+00 0.00037 ];
INF_CAPT                  (idx, [1:   4]) = [  7.55849E-03 0.00087  2.91979E-02 0.00060 ];
INF_ABS                   (idx, [1:   4]) = [  1.14236E-02 0.00071  1.14159E-01 0.00077 ];
INF_FISS                  (idx, [1:   4]) = [  3.86506E-03 0.00101  8.49613E-02 0.00083 ];
INF_NSF                   (idx, [1:   4]) = [  9.75630E-03 0.00099  2.06983E-01 0.00083 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52423E+00 8.9E-05  2.43620E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03060E+02 7.2E-06  2.02270E+02 2.7E-09 ];
INF_INVV                  (idx, [1:   4]) = [  5.49785E-08 0.00094  2.27396E-06 0.00042 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34752E-01 0.00019  1.22301E+00 0.00043 ];
INF_SCATT1                (idx, [1:   4]) = [  2.36126E-01 0.00027  3.30454E-01 0.00116 ];
INF_SCATT2                (idx, [1:   4]) = [  9.35508E-02 0.00048  8.41861E-02 0.00309 ];
INF_SCATT3                (idx, [1:   4]) = [  7.30558E-03 0.00525  2.56251E-02 0.00797 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.55485E-03 0.00472 -5.15604E-03 0.02739 ];
INF_SCATT5                (idx, [1:   4]) = [  4.44978E-04 0.08795  4.57992E-03 0.03254 ];
INF_SCATT6                (idx, [1:   4]) = [  4.99009E-03 0.00627 -1.23045E-02 0.01259 ];
INF_SCATT7                (idx, [1:   4]) = [  7.28347E-04 0.04360 -5.85369E-04 0.18594 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34792E-01 0.00019  1.22301E+00 0.00043 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.36126E-01 0.00027  3.30454E-01 0.00116 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.35511E-02 0.00048  8.41861E-02 0.00309 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.30536E-03 0.00527  2.56251E-02 0.00797 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.55501E-03 0.00472 -5.15604E-03 0.02739 ];
INF_SCATTP5               (idx, [1:   4]) = [  4.44839E-04 0.08794  4.57992E-03 0.03254 ];
INF_SCATTP6               (idx, [1:   4]) = [  4.98964E-03 0.00627 -1.23045E-02 0.01259 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.28371E-04 0.04366 -5.85369E-04 0.18594 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30844E-01 0.00050  8.85495E-01 0.00044 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44398E+00 0.00050  3.76439E-01 0.00044 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13828E-02 0.00071  1.14159E-01 0.00077 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72498E-02 0.00038  1.16520E-01 0.00083 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18922E-01 0.00019  1.58301E-02 0.00079  2.26020E-03 0.00796  1.22075E+00 0.00043 ];
INF_S1                    (idx, [1:   8]) = [  2.31552E-01 0.00028  4.57434E-03 0.00177  9.09607E-04 0.01818  3.29544E-01 0.00115 ];
INF_S2                    (idx, [1:   8]) = [  9.49276E-02 0.00048 -1.37682E-03 0.00534  4.87970E-04 0.02420  8.36981E-02 0.00307 ];
INF_S3                    (idx, [1:   8]) = [  8.92372E-03 0.00422 -1.61814E-03 0.00356  1.68901E-04 0.06406  2.54562E-02 0.00782 ];
INF_S4                    (idx, [1:   8]) = [ -9.03345E-03 0.00512 -5.21397E-04 0.01001 -7.54368E-06 0.95336 -5.14849E-03 0.02752 ];
INF_S5                    (idx, [1:   8]) = [  4.17455E-04 0.09721  2.75235E-05 0.24606 -7.76213E-05 0.08770  4.65754E-03 0.03190 ];
INF_S6                    (idx, [1:   8]) = [  5.11622E-03 0.00606 -1.26131E-04 0.02331 -9.82431E-05 0.06672 -1.22062E-02 0.01284 ];
INF_S7                    (idx, [1:   8]) = [  8.82798E-04 0.03487 -1.54451E-04 0.02483 -8.27719E-05 0.06294 -5.02597E-04 0.21665 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18962E-01 0.00019  1.58301E-02 0.00079  2.26020E-03 0.00796  1.22075E+00 0.00043 ];
INF_SP1                   (idx, [1:   8]) = [  2.31552E-01 0.00028  4.57434E-03 0.00177  9.09607E-04 0.01818  3.29544E-01 0.00115 ];
INF_SP2                   (idx, [1:   8]) = [  9.49280E-02 0.00048 -1.37682E-03 0.00534  4.87970E-04 0.02420  8.36981E-02 0.00307 ];
INF_SP3                   (idx, [1:   8]) = [  8.92351E-03 0.00423 -1.61814E-03 0.00356  1.68901E-04 0.06406  2.54562E-02 0.00782 ];
INF_SP4                   (idx, [1:   8]) = [ -9.03362E-03 0.00512 -5.21397E-04 0.01001 -7.54368E-06 0.95336 -5.14849E-03 0.02752 ];
INF_SP5                   (idx, [1:   8]) = [  4.17316E-04 0.09722  2.75235E-05 0.24606 -7.76213E-05 0.08770  4.65754E-03 0.03190 ];
INF_SP6                   (idx, [1:   8]) = [  5.11577E-03 0.00606 -1.26131E-04 0.02331 -9.82431E-05 0.06672 -1.22062E-02 0.01284 ];
INF_SP7                   (idx, [1:   8]) = [  8.82821E-04 0.03492 -1.54451E-04 0.02483 -8.27719E-05 0.06294 -5.02597E-04 0.21665 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42518E-01 0.00125  7.95535E-01 0.00652 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.43055E-01 0.00189  8.04263E-01 0.01101 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.42485E-01 0.00209  7.88135E-01 0.01313 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.42058E-01 0.00199  7.99414E-01 0.00960 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37452E+00 0.00124  4.19440E-01 0.00663 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37155E+00 0.00189  4.15728E-01 0.01160 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37480E+00 0.00209  4.24669E-01 0.01297 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37721E+00 0.00199  4.17922E-01 0.00989 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  6.99831E-03 0.01642  2.04040E-04 0.10772  9.40448E-04 0.04413  7.03475E-04 0.05497  1.36377E-03 0.03634  2.22262E-03 0.02795  7.08725E-04 0.05045  6.05005E-04 0.05502  2.50219E-04 0.09532 ];
LAMBDA                    (idx, [1:  18]) = [  4.56703E-01 0.02656  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.0E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Wed Feb  7 19:56:16 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 9.3E-10  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.68067E-02 0.00143  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73193E-01 3.9E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.33941E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98749E-01 6.9E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98240E-01 9.6E-06 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.33617E-01 0.00013  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.83434E+00 0.00049  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.70085E+01 0.00053  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.70085E+01 0.00053  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.76122E+00 0.00069  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.26023E-01 0.00157  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000465 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00093E+03 0.00126 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00093E+03 0.00126 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  6.34310E+01 ;
RUNNING_TIME              (idx, 1)        =  6.68677E+01 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  3.69367E-01  2.01367E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  6.45773E+01  2.12083E+01  2.37394E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  5.87333E-02  5.35333E-02 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  8.18333E-03  9.00000E-04 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  6.68672E+01  2.43632E+03 ];
CPU_USAGE                 (idx, 1)        = 0.94860 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  1.00005E+00 0.00010 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.25666E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  2.40000E-05 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 59 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  3.17467E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.23306E+01 ;
TOT_SF_RATE               (idx, 1)        =  4.05571E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  3.13937E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  2.30223E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  2.86073E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.21003E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  1.39466E+03 ;
INGESTION_TOXICITY        (idx, 1)        =  2.55659E+03 ;
ACTINIDE_INH_TOX          (idx, 1)        =  3.02932E+02 ;
ACTINIDE_ING_TOX          (idx, 1)        =  2.58306E+02 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  1.09172E+03 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  2.29828E+03 ;
SR90_ACTIVITY             (idx, 1)        =  7.03808E+06 ;
TE132_ACTIVITY            (idx, 1)        =  1.73415E+10 ;
I131_ACTIVITY             (idx, 1)        =  3.52286E+09 ;
I132_ACTIVITY             (idx, 1)        =  1.13533E+10 ;
CS134_ACTIVITY            (idx, 1)        =  3.92394E+02 ;
CS137_ACTIVITY            (idx, 1)        =  7.48171E+06 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  3.39771E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.96871E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.33798E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  4.11076E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.71745E+09 0.00100  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 1 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  1.00000E-02  1.00036E-02 ];
BURN_DAYS                 (idx, [1:  2])  = [  2.97885E-01  2.97885E-01 ];
FIMA                      (idx, [1:  3])  = [  1.05259E-05  1.67283E+17  1.58924E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.86260E-01 0.00236 ];
U235_FISS                 (idx, [1:   4]) = [  6.15002E+12 0.00106  9.46263E-01 0.00033 ];
U238_FISS                 (idx, [1:   4]) = [  3.48445E+11 0.00612  5.35817E-02 0.00584 ];
PU239_FISS                (idx, [1:   4]) = [  5.61869E+07 0.44567  8.63134E-06 0.44578 ];
U235_CAPT                 (idx, [1:   4]) = [  1.46036E+12 0.00274  2.94945E-01 0.00234 ];
U238_CAPT                 (idx, [1:   4]) = [  2.90269E+12 0.00229  5.86002E-01 0.00124 ];
PU239_CAPT                (idx, [1:   4]) = [  3.52702E+07 0.57653  7.01502E-06 0.57642 ];
XE135_CAPT                (idx, [1:   4]) = [  7.42255E+10 0.01283  1.49922E-02 0.01275 ];
SM149_CAPT                (idx, [1:   4]) = [  6.83047E+07 0.40641  1.36075E-05 0.40634 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000465 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51955E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000465 1.00152E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 432521 4.33037E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 567944 5.68482E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000465 1.00152E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -1.51340E-09 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.59734E+13 2.2E-05  1.59734E+13 2.2E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49722E+12 1.9E-06  6.49722E+12 1.9E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  4.95100E+12 0.00097  4.49506E+12 0.00106  4.55938E+11 0.00119 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.14482E+13 0.00042  1.09923E+13 0.00043  4.55938E+11 0.00119 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.14349E+13 0.00100  1.14349E+13 0.00100  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.82604E+14 0.00085  1.60272E+14 0.00093  3.22332E+14 0.00088 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.14482E+13 0.00042 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.08849E+14 0.00071 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27830E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27830E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.88115E+00 0.00074 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.52595E-01 0.00030 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.76876E-01 0.00085 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.35241E+00 0.00080 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.39763E+00 0.00088 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.39763E+00 0.00088 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45850E+00 2.4E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02470E+02 1.9E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.39716E+00 0.00093  1.38790E+00 0.00089  9.72316E-03 0.01668 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.39751E+00 0.00042 ];
COL_KEFF                  (idx, [1:   2]) = [  1.39759E+00 0.00100 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.39751E+00 0.00042 ];
ABS_KINF                  (idx, [1:   2]) = [  1.39751E+00 0.00042 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.70295E+01 0.00039 ];
IMP_ALF                   (idx, [1:   2]) = [  1.70362E+01 0.00018 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.13046E-07 0.00681 ];
IMP_EALF                  (idx, [1:   2]) = [  8.00353E-07 0.00305 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.90976E-01 0.00611 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.90627E-01 0.00238 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.23546E-03 0.01204  1.44135E-04 0.07217  7.76817E-04 0.03023  4.70466E-04 0.04054  9.99827E-04 0.02714  1.66350E-03 0.02132  5.56695E-04 0.03628  4.70411E-04 0.03915  1.53605E-04 0.07116 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.50315E-01 0.02008  4.03921E-03 0.06466  2.52928E-02 0.01541  3.09578E-02 0.02736  1.23197E-01 0.01265  2.87203E-01 0.00606  5.23859E-01 0.02336  1.19012E+00 0.02736  1.20145E+00 0.06265 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.21116E-03 0.01763  1.87946E-04 0.10381  1.05052E-03 0.04289  6.80467E-04 0.05404  1.34797E-03 0.03863  2.33053E-03 0.03273  7.55580E-04 0.05129  6.55911E-04 0.05736  2.02239E-04 0.09194 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.41606E-01 0.02538  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.7E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.22600E-05 0.00207  1.22512E-05 0.00208  1.35090E-05 0.02177 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.71217E-05 0.00185  1.71094E-05 0.00185  1.88677E-05 0.02177 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.95128E-03 0.01695  1.86307E-04 0.10331  1.02588E-03 0.04245  6.53151E-04 0.05460  1.27695E-03 0.03803  2.23574E-03 0.02984  7.05093E-04 0.05067  6.46259E-04 0.05497  2.21901E-04 0.09408 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.52690E-01 0.02761  1.24667E-02 0.0E+00  2.82917E-02 1.1E-09  4.25244E-02 8.0E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.5E-09  3.55460E+00 5.1E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.22439E-05 0.00438  1.22316E-05 0.00440  1.04383E-05 0.05213 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.70998E-05 0.00430  1.70828E-05 0.00432  1.45805E-05 0.05194 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  6.97965E-03 0.04844  1.89637E-04 0.33673  1.08866E-03 0.11943  7.46878E-04 0.15755  1.18516E-03 0.10770  2.04991E-03 0.08877  8.31362E-04 0.15658  6.06558E-04 0.17178  2.81480E-04 0.26952 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.80230E-01 0.07014  1.24667E-02 5.5E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.0E-09  2.92467E-01 5.8E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 6.0E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.06936E-03 0.04652  2.18291E-04 0.31863  1.12660E-03 0.11909  7.13678E-04 0.14932  1.21316E-03 0.10343  2.01810E-03 0.08519  8.39014E-04 0.14503  6.15481E-04 0.16553  3.25030E-04 0.25495 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.92496E-01 0.06972  1.24667E-02 3.9E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.0E-09  2.92467E-01 5.8E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 6.0E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -5.74524E+02 0.04893 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.22845E-05 0.00124 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.71560E-05 0.00083 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.74182E-03 0.00833 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.49421E+02 0.00854 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.18179E-07 0.00117 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92793E-06 0.00108  2.92787E-06 0.00109  2.92900E-06 0.01233 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.97329E-05 0.00129  1.97294E-05 0.00130  2.01377E-05 0.01547 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.77504E-01 0.00085  5.76120E-01 0.00086  9.01772E-01 0.01956 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.27439E+01 0.02718 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.69747E+01 0.00053  2.96213E+01 0.00066 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.30891E+04 0.00765  5.43614E+04 0.00292  1.14134E+05 0.00252  1.25749E+05 0.00173  1.18468E+05 0.00157  1.32834E+05 0.00167  9.07340E+04 0.00137  8.12204E+04 0.00132  6.18069E+04 0.00174  5.02202E+04 0.00151  4.31460E+04 0.00172  3.90757E+04 0.00187  3.59498E+04 0.00213  3.41254E+04 0.00217  3.30753E+04 0.00187  2.85713E+04 0.00235  2.80550E+04 0.00221  2.78402E+04 0.00205  2.72300E+04 0.00245  5.27466E+04 0.00128  5.01886E+04 0.00160  3.57102E+04 0.00179  2.27895E+04 0.00284  2.58200E+04 0.00262  2.39049E+04 0.00253  2.17172E+04 0.00247  3.48552E+04 0.00170  8.02058E+03 0.00446  1.01001E+04 0.00360  9.19157E+03 0.00368  5.30314E+03 0.00359  9.31311E+03 0.00418  6.34675E+03 0.00414  5.27343E+03 0.00463  9.95135E+02 0.00908  9.90299E+02 0.00868  1.00032E+03 0.01041  1.03296E+03 0.00840  1.03820E+03 0.00669  1.01400E+03 0.00814  1.05254E+03 0.00878  9.84885E+02 0.00968  1.88669E+03 0.00556  3.01497E+03 0.00662  3.82758E+03 0.00392  1.00960E+04 0.00360  1.04303E+04 0.00249  1.08663E+04 0.00400  6.92662E+03 0.00389  4.82395E+03 0.00478  3.63840E+03 0.00525  4.07750E+03 0.00420  7.29742E+03 0.00402  9.04347E+03 0.00295  1.63388E+04 0.00233  2.22374E+04 0.00257  2.97791E+04 0.00186  1.76857E+04 0.00234  1.20774E+04 0.00248  8.41387E+03 0.00310  7.35687E+03 0.00358  7.03086E+03 0.00331  5.80581E+03 0.00339  3.84681E+03 0.00320  3.52111E+03 0.00508  3.04866E+03 0.00458  2.54995E+03 0.00462  1.97179E+03 0.00437  1.28315E+03 0.00729  4.36906E+02 0.00859 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.39803E+00 0.00100 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.25349E+14 0.00096  5.73128E+13 0.00069 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.45872E-01 0.00033  1.33865E+00 0.00060 ];
INF_CAPT                  (idx, [1:   4]) = [  7.54209E-03 0.00119  3.04276E-02 0.00051 ];
INF_ABS                   (idx, [1:   4]) = [  1.14066E-02 0.00097  1.15167E-01 0.00063 ];
INF_FISS                  (idx, [1:   4]) = [  3.86450E-03 0.00099  8.47395E-02 0.00067 ];
INF_NSF                   (idx, [1:   4]) = [  9.75525E-03 0.00097  2.06443E-01 0.00067 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52432E+00 7.9E-05  2.43620E+00 0.0E+00 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03060E+02 7.1E-06  2.02270E+02 0.0E+00 ];
INF_INVV                  (idx, [1:   4]) = [  5.50196E-08 0.00106  2.27132E-06 0.00052 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34457E-01 0.00034  1.22361E+00 0.00069 ];
INF_SCATT1                (idx, [1:   4]) = [  2.35927E-01 0.00046  3.30813E-01 0.00107 ];
INF_SCATT2                (idx, [1:   4]) = [  9.35542E-02 0.00053  8.46662E-02 0.00295 ];
INF_SCATT3                (idx, [1:   4]) = [  7.41063E-03 0.00521  2.58982E-02 0.00806 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.35916E-03 0.00493 -5.33274E-03 0.03151 ];
INF_SCATT5                (idx, [1:   4]) = [  5.16690E-04 0.07663  4.46267E-03 0.03884 ];
INF_SCATT6                (idx, [1:   4]) = [  5.02957E-03 0.00734 -1.23020E-02 0.01047 ];
INF_SCATT7                (idx, [1:   4]) = [  7.22137E-04 0.04399 -6.36477E-04 0.20151 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34498E-01 0.00034  1.22361E+00 0.00069 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.35926E-01 0.00046  3.30813E-01 0.00107 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.35541E-02 0.00053  8.46662E-02 0.00295 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.41001E-03 0.00521  2.58982E-02 0.00806 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.35964E-03 0.00493 -5.33274E-03 0.03151 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.16497E-04 0.07653  4.46267E-03 0.03884 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.02954E-03 0.00731 -1.23020E-02 0.01047 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.22305E-04 0.04403 -6.36477E-04 0.20151 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30790E-01 0.00067  8.86334E-01 0.00075 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44433E+00 0.00067  3.76086E-01 0.00076 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13657E-02 0.00098  1.15167E-01 0.00063 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72222E-02 0.00051  1.17253E-01 0.00118 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18650E-01 0.00033  1.58068E-02 0.00096  2.21352E-03 0.01126  1.22140E+00 0.00070 ];
INF_S1                    (idx, [1:   8]) = [  2.31377E-01 0.00046  4.55018E-03 0.00215  8.88665E-04 0.01698  3.29924E-01 0.00106 ];
INF_S2                    (idx, [1:   8]) = [  9.49557E-02 0.00050 -1.40154E-03 0.00594  4.56803E-04 0.02384  8.42094E-02 0.00298 ];
INF_S3                    (idx, [1:   8]) = [  9.02948E-03 0.00397 -1.61885E-03 0.00424  1.55545E-04 0.05146  2.57426E-02 0.00810 ];
INF_S4                    (idx, [1:   8]) = [ -8.83393E-03 0.00491 -5.25229E-04 0.01177 -1.28567E-05 0.50949 -5.31988E-03 0.03158 ];
INF_S5                    (idx, [1:   8]) = [  4.91678E-04 0.07868  2.50115E-05 0.20942 -7.70941E-05 0.07038  4.53976E-03 0.03866 ];
INF_S6                    (idx, [1:   8]) = [  5.15016E-03 0.00698 -1.20592E-04 0.03803 -9.23656E-05 0.06713 -1.22097E-02 0.01049 ];
INF_S7                    (idx, [1:   8]) = [  8.71642E-04 0.03713 -1.49504E-04 0.02653 -7.13266E-05 0.07476 -5.65150E-04 0.22655 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18691E-01 0.00033  1.58068E-02 0.00096  2.21352E-03 0.01126  1.22140E+00 0.00070 ];
INF_SP1                   (idx, [1:   8]) = [  2.31376E-01 0.00046  4.55018E-03 0.00215  8.88665E-04 0.01698  3.29924E-01 0.00106 ];
INF_SP2                   (idx, [1:   8]) = [  9.49556E-02 0.00050 -1.40154E-03 0.00594  4.56803E-04 0.02384  8.42094E-02 0.00298 ];
INF_SP3                   (idx, [1:   8]) = [  9.02886E-03 0.00397 -1.61885E-03 0.00424  1.55545E-04 0.05146  2.57426E-02 0.00810 ];
INF_SP4                   (idx, [1:   8]) = [ -8.83441E-03 0.00491 -5.25229E-04 0.01177 -1.28567E-05 0.50949 -5.31988E-03 0.03158 ];
INF_SP5                   (idx, [1:   8]) = [  4.91486E-04 0.07855  2.50115E-05 0.20942 -7.70941E-05 0.07038  4.53976E-03 0.03866 ];
INF_SP6                   (idx, [1:   8]) = [  5.15013E-03 0.00695 -1.20592E-04 0.03803 -9.23656E-05 0.06713 -1.22097E-02 0.01049 ];
INF_SP7                   (idx, [1:   8]) = [  8.71810E-04 0.03717 -1.49504E-04 0.02653 -7.13266E-05 0.07476 -5.65150E-04 0.22655 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42706E-01 0.00103  7.97545E-01 0.00620 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.42873E-01 0.00206  7.99503E-01 0.01155 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.42606E-01 0.00168  7.98051E-01 0.01041 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.42680E-01 0.00170  7.99252E-01 0.01000 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37344E+00 0.00104  4.18333E-01 0.00617 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37260E+00 0.00205  4.18218E-01 0.01121 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37407E+00 0.00169  4.18730E-01 0.01004 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37365E+00 0.00170  4.18051E-01 0.00993 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.21116E-03 0.01763  1.87946E-04 0.10381  1.05052E-03 0.04289  6.80467E-04 0.05404  1.34797E-03 0.03863  2.33053E-03 0.03273  7.55580E-04 0.05129  6.55911E-04 0.05736  2.02239E-04 0.09194 ];
LAMBDA                    (idx, [1:  18]) = [  4.41606E-01 0.02538  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.7E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Wed Feb  7 20:43:53 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 2.0E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.68561E-02 0.00136  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73144E-01 3.8E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.34518E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98766E-01 7.5E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98262E-01 1.1E-05 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.34205E-01 0.00013  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.82964E+00 0.00050  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.68831E+01 0.00050  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.68831E+01 0.00050  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.68725E+00 0.00070  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.22588E-01 0.00151  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000597 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00119E+03 0.00140 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00119E+03 0.00140 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.10683E+02 ;
RUNNING_TIME              (idx, 1)        =  1.14490E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  5.96300E-01  1.24083E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  1.11960E+02  2.68574E+01  2.05253E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  6.84833E-02  4.80000E-03 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  1.63500E-02  9.16665E-04 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  1.14489E+02  2.34003E+03 ];
CPU_USAGE                 (idx, 1)        = 0.96676 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.99452E-01 0.00010 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.53323E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  2.40000E-05 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 78 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  3.51189E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.28627E+01 ;
TOT_SF_RATE               (idx, 1)        =  4.05689E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  3.79480E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  2.77090E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  3.13240E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.25856E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  3.47658E+03 ;
INGESTION_TOXICITY        (idx, 1)        =  5.74319E+03 ;
ACTINIDE_INH_TOX          (idx, 1)        =  9.35205E+02 ;
ACTINIDE_ING_TOX          (idx, 1)        =  7.59592E+02 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  2.54137E+03 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  4.98360E+03 ;
SR90_ACTIVITY             (idx, 1)        =  2.83960E+07 ;
TE132_ACTIVITY            (idx, 1)        =  6.34899E+10 ;
I131_ACTIVITY             (idx, 1)        =  1.59937E+10 ;
I132_ACTIVITY             (idx, 1)        =  5.81332E+10 ;
CS134_ACTIVITY            (idx, 1)        =  2.48760E+03 ;
CS137_ACTIVITY            (idx, 1)        =  3.02328E+07 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  3.73652E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.97276E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.34052E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  4.79353E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.83390E+09 0.00097  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 2 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  4.00000E-02  4.00150E-02 ];
BURN_DAYS                 (idx, [1:  2])  = [  1.19154E+00  8.93655E-01 ];
FIMA                      (idx, [1:  3])  = [  4.21052E-05  6.69160E+17  1.58919E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.92373E-01 0.00229 ];
U235_FISS                 (idx, [1:   4]) = [  6.14199E+12 0.00115  9.45200E-01 0.00032 ];
U238_FISS                 (idx, [1:   4]) = [  3.54387E+11 0.00580  5.45077E-02 0.00550 ];
PU239_FISS                (idx, [1:   4]) = [  8.85952E+08 0.11050  1.36283E-04 0.11044 ];
U235_CAPT                 (idx, [1:   4]) = [  1.48080E+12 0.00294  2.85432E-01 0.00236 ];
U238_CAPT                 (idx, [1:   4]) = [  2.95526E+12 0.00216  5.69634E-01 0.00126 ];
PU239_CAPT                (idx, [1:   4]) = [  7.18683E+08 0.12921  1.38099E-04 0.12945 ];
XE135_CAPT                (idx, [1:   4]) = [  2.33068E+11 0.00652  4.49763E-02 0.00664 ];
SM149_CAPT                (idx, [1:   4]) = [  1.20214E+09 0.09509  2.31648E-04 0.09504 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000597 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.56405E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000597 1.00156E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 444102 4.44558E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 556495 5.57006E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000597 1.00156E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -3.49246E-10 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.59770E+13 2.3E-05  1.59770E+13 2.3E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49705E+12 2.0E-06  6.49705E+12 2.0E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  5.17466E+12 0.00094  4.71549E+12 0.00101  4.59172E+11 0.00121 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.16717E+13 0.00042  1.12125E+13 0.00043  4.59172E+11 0.00121 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.16678E+13 0.00097  1.16678E+13 0.00097  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.91309E+14 0.00083  1.63251E+14 0.00094  3.28058E+14 0.00084 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.16717E+13 0.00042 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.13699E+14 0.00070 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27811E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27811E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.83039E+00 0.00077 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.53970E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.76672E-01 0.00085 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36065E+00 0.00079 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.36972E+00 0.00094 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.36972E+00 0.00094 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45912E+00 2.5E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02475E+02 2.0E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.36999E+00 0.00100  1.36005E+00 0.00095  9.67012E-03 0.01603 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.37114E+00 0.00041 ];
COL_KEFF                  (idx, [1:   2]) = [  1.36996E+00 0.00096 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.37114E+00 0.00041 ];
ABS_KINF                  (idx, [1:   2]) = [  1.37114E+00 0.00041 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.69820E+01 0.00036 ];
IMP_ALF                   (idx, [1:   2]) = [  1.69825E+01 0.00018 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.51160E-07 0.00626 ];
IMP_EALF                  (idx, [1:   2]) = [  8.44601E-07 0.00311 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.94905E-01 0.00577 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.95239E-01 0.00242 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.27370E-03 0.01186  1.48830E-04 0.07009  7.62512E-04 0.02923  4.39481E-04 0.03913  1.00819E-03 0.02732  1.73116E-03 0.01954  5.52709E-04 0.03783  4.80264E-04 0.03972  1.50556E-04 0.07086 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.48435E-01 0.01935  4.18881E-03 0.06293  2.58586E-02 0.01373  3.01923E-02 0.02861  1.25326E-01 0.01111  2.91297E-01 0.00284  5.23859E-01 0.02336  1.16723E+00 0.02833  1.16591E+00 0.06408 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.14232E-03 0.01606  2.01845E-04 0.09367  1.03902E-03 0.04194  5.91286E-04 0.05295  1.34812E-03 0.03750  2.36668E-03 0.02649  7.43334E-04 0.05486  6.34224E-04 0.05345  2.17817E-04 0.09921 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.51888E-01 0.02672  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 6.0E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.0E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.21592E-05 0.00203  1.21540E-05 0.00204  1.29211E-05 0.02197 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.66491E-05 0.00172  1.66421E-05 0.00174  1.76847E-05 0.02184 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.06918E-03 0.01632  1.90071E-04 0.10447  1.01836E-03 0.03855  6.05188E-04 0.05559  1.30352E-03 0.03862  2.37386E-03 0.02824  7.28506E-04 0.05473  6.47663E-04 0.05273  2.02009E-04 0.09504 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.59039E-01 0.02768  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 8.1E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.5E-09  3.55460E+00 5.1E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.22402E-05 0.00447  1.22402E-05 0.00449  9.39730E-06 0.04946 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.67635E-05 0.00444  1.67636E-05 0.00447  1.28337E-05 0.04933 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.40439E-03 0.04845  3.08693E-04 0.23791  1.01016E-03 0.12454  5.58898E-04 0.17979  1.43579E-03 0.10801  2.33001E-03 0.08503  8.27983E-04 0.14747  6.79400E-04 0.15204  2.53454E-04 0.24859 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.81721E-01 0.06570  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.4E-09  2.92467E-01 6.0E-09  6.66488E-01 5.0E-09  1.63478E+00 0.0E+00  3.55460E+00 5.4E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.38625E-03 0.04655  3.20818E-04 0.22260  9.98036E-04 0.11943  5.56694E-04 0.17526  1.42899E-03 0.10633  2.34569E-03 0.08294  8.14735E-04 0.14178  6.72289E-04 0.14583  2.48992E-04 0.23761 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.76567E-01 0.06513  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.2E-09  2.92467E-01 5.9E-09  6.66488E-01 5.0E-09  1.63478E+00 0.0E+00  3.55460E+00 6.0E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -6.14067E+02 0.04902 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.21978E-05 0.00132 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.67023E-05 0.00081 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.10959E-03 0.01033 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.83189E+02 0.01037 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.12658E-07 0.00112 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92253E-06 0.00101  2.92257E-06 0.00101  2.91543E-06 0.01303 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.92962E-05 0.00122  1.92969E-05 0.00122  1.92882E-05 0.01561 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.77365E-01 0.00086  5.76000E-01 0.00087  9.01927E-01 0.02157 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.20571E+01 0.02766 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.68499E+01 0.00050  2.93595E+01 0.00062 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.33263E+04 0.00777  5.40794E+04 0.00335  1.14186E+05 0.00165  1.25569E+05 0.00168  1.18474E+05 0.00164  1.33268E+05 0.00174  9.04076E+04 0.00199  8.14274E+04 0.00147  6.17373E+04 0.00162  5.03693E+04 0.00233  4.31344E+04 0.00181  3.89183E+04 0.00212  3.58947E+04 0.00191  3.41635E+04 0.00254  3.31770E+04 0.00218  2.85445E+04 0.00197  2.81155E+04 0.00235  2.78689E+04 0.00260  2.72104E+04 0.00212  5.26547E+04 0.00137  5.02336E+04 0.00208  3.58163E+04 0.00135  2.27560E+04 0.00264  2.57641E+04 0.00201  2.39875E+04 0.00216  2.16948E+04 0.00281  3.47779E+04 0.00150  8.01299E+03 0.00395  1.01612E+04 0.00361  9.23796E+03 0.00441  5.29121E+03 0.00381  9.33736E+03 0.00418  6.27601E+03 0.00422  5.27451E+03 0.00396  9.73896E+02 0.00791  9.69209E+02 0.00809  9.99075E+02 0.00875  1.03343E+03 0.00773  1.02673E+03 0.00745  1.01512E+03 0.00813  1.04839E+03 0.00694  9.86126E+02 0.00720  1.86286E+03 0.00569  2.96932E+03 0.00515  3.80454E+03 0.00497  1.00056E+04 0.00348  1.04017E+04 0.00322  1.08926E+04 0.00342  6.94411E+03 0.00348  4.75783E+03 0.00358  3.57580E+03 0.00495  4.05760E+03 0.00428  7.25994E+03 0.00389  8.95483E+03 0.00363  1.60542E+04 0.00256  2.17761E+04 0.00251  2.91553E+04 0.00252  1.71607E+04 0.00252  1.17175E+04 0.00242  8.16773E+03 0.00279  7.14657E+03 0.00304  6.85648E+03 0.00304  5.61900E+03 0.00389  3.72012E+03 0.00460  3.40822E+03 0.00439  2.97716E+03 0.00494  2.48140E+03 0.00440  1.91998E+03 0.00490  1.23887E+03 0.00637  4.28089E+02 0.00899 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.37017E+00 0.00064 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.33996E+14 0.00070  5.73673E+13 0.00096 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.45880E-01 0.00031  1.33949E+00 0.00039 ];
INF_CAPT                  (idx, [1:   4]) = [  7.55223E-03 0.00147  3.30842E-02 0.00074 ];
INF_ABS                   (idx, [1:   4]) = [  1.14134E-02 0.00105  1.17181E-01 0.00089 ];
INF_FISS                  (idx, [1:   4]) = [  3.86113E-03 0.00081  8.40968E-02 0.00095 ];
INF_NSF                   (idx, [1:   4]) = [  9.74877E-03 0.00080  2.04883E-01 0.00095 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52485E+00 7.5E-05  2.43628E+00 3.9E-08 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03062E+02 5.4E-06  2.02271E+02 3.8E-09 ];
INF_INVV                  (idx, [1:   4]) = [  5.48913E-08 0.00081  2.26279E-06 0.00034 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34456E-01 0.00032  1.22224E+00 0.00044 ];
INF_SCATT1                (idx, [1:   4]) = [  2.35992E-01 0.00045  3.31214E-01 0.00108 ];
INF_SCATT2                (idx, [1:   4]) = [  9.35425E-02 0.00071  8.47235E-02 0.00351 ];
INF_SCATT3                (idx, [1:   4]) = [  7.36202E-03 0.00719  2.57845E-02 0.00690 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.40799E-03 0.00468 -5.19547E-03 0.02993 ];
INF_SCATT5                (idx, [1:   4]) = [  5.42746E-04 0.09644  4.61991E-03 0.03158 ];
INF_SCATT6                (idx, [1:   4]) = [  5.12275E-03 0.00662 -1.22606E-02 0.01311 ];
INF_SCATT7                (idx, [1:   4]) = [  7.75276E-04 0.04718 -4.44760E-04 0.29354 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34498E-01 0.00032  1.22224E+00 0.00044 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.35993E-01 0.00045  3.31214E-01 0.00108 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.35424E-02 0.00071  8.47235E-02 0.00351 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.36226E-03 0.00721  2.57845E-02 0.00690 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.40733E-03 0.00468 -5.19547E-03 0.02993 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.42755E-04 0.09648  4.61991E-03 0.03158 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.12323E-03 0.00662 -1.22606E-02 0.01311 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.74923E-04 0.04716 -4.44760E-04 0.29354 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30833E-01 0.00063  8.85976E-01 0.00044 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44406E+00 0.00063  3.76234E-01 0.00044 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13713E-02 0.00107  1.17181E-01 0.00089 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72343E-02 0.00047  1.19565E-01 0.00095 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18646E-01 0.00031  1.58100E-02 0.00097  2.31710E-03 0.00913  1.21992E+00 0.00045 ];
INF_S1                    (idx, [1:   8]) = [  2.31430E-01 0.00045  4.56274E-03 0.00230  9.30246E-04 0.01422  3.30284E-01 0.00108 ];
INF_S2                    (idx, [1:   8]) = [  9.49142E-02 0.00071 -1.37173E-03 0.00599  5.07224E-04 0.02614  8.42163E-02 0.00354 ];
INF_S3                    (idx, [1:   8]) = [  8.97413E-03 0.00606 -1.61211E-03 0.00459  1.76450E-04 0.05190  2.56080E-02 0.00697 ];
INF_S4                    (idx, [1:   8]) = [ -8.87784E-03 0.00508 -5.30150E-04 0.01093 -1.00344E-05 0.81586 -5.18544E-03 0.03034 ];
INF_S5                    (idx, [1:   8]) = [  5.21857E-04 0.09836  2.08892E-05 0.26908 -7.14432E-05 0.10017  4.69136E-03 0.03078 ];
INF_S6                    (idx, [1:   8]) = [  5.24836E-03 0.00633 -1.25612E-04 0.04173 -8.99714E-05 0.06064 -1.21707E-02 0.01310 ];
INF_S7                    (idx, [1:   8]) = [  9.32377E-04 0.03870 -1.57101E-04 0.02854 -7.69972E-05 0.05778 -3.67763E-04 0.35849 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18688E-01 0.00031  1.58100E-02 0.00097  2.31710E-03 0.00913  1.21992E+00 0.00045 ];
INF_SP1                   (idx, [1:   8]) = [  2.31430E-01 0.00045  4.56274E-03 0.00230  9.30246E-04 0.01422  3.30284E-01 0.00108 ];
INF_SP2                   (idx, [1:   8]) = [  9.49142E-02 0.00071 -1.37173E-03 0.00599  5.07224E-04 0.02614  8.42163E-02 0.00354 ];
INF_SP3                   (idx, [1:   8]) = [  8.97437E-03 0.00607 -1.61211E-03 0.00459  1.76450E-04 0.05190  2.56080E-02 0.00697 ];
INF_SP4                   (idx, [1:   8]) = [ -8.87718E-03 0.00508 -5.30150E-04 0.01093 -1.00344E-05 0.81586 -5.18544E-03 0.03034 ];
INF_SP5                   (idx, [1:   8]) = [  5.21866E-04 0.09840  2.08892E-05 0.26908 -7.14432E-05 0.10017  4.69136E-03 0.03078 ];
INF_SP6                   (idx, [1:   8]) = [  5.24884E-03 0.00633 -1.25612E-04 0.04173 -8.99714E-05 0.06064 -1.21707E-02 0.01310 ];
INF_SP7                   (idx, [1:   8]) = [  9.32024E-04 0.03869 -1.57101E-04 0.02854 -7.69972E-05 0.05778 -3.67763E-04 0.35849 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42779E-01 0.00130  7.94316E-01 0.00524 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.43492E-01 0.00182  8.01062E-01 0.00849 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.42923E-01 0.00229  7.96319E-01 0.01001 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.41974E-01 0.00212  7.88427E-01 0.00718 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37305E+00 0.00129  4.19922E-01 0.00518 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.36908E+00 0.00181  4.16845E-01 0.00861 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37235E+00 0.00228  4.19609E-01 0.01012 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37771E+00 0.00213  4.23311E-01 0.00725 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.14232E-03 0.01606  2.01845E-04 0.09367  1.03902E-03 0.04194  5.91286E-04 0.05295  1.34812E-03 0.03750  2.36668E-03 0.02649  7.43334E-04 0.05486  6.34224E-04 0.05345  2.17817E-04 0.09921 ];
LAMBDA                    (idx, [1:  18]) = [  4.51888E-01 0.02672  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 6.0E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.0E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Wed Feb  7 21:43:01 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.69202E-02 0.00140  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73080E-01 3.9E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.34419E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98750E-01 7.0E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98240E-01 9.9E-06 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.34092E-01 0.00013  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.82793E+00 0.00050  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.68764E+01 0.00057  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.68764E+01 0.00057  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.68973E+00 0.00070  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.25261E-01 0.00151  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000384 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00077E+03 0.00131 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00077E+03 0.00131 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  1.69634E+02 ;
RUNNING_TIME              (idx, 1)        =  1.73627E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  9.57833E-01  1.91333E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  1.70723E+02  2.92029E+01  2.95600E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  8.00167E-02  5.78333E-03 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  2.62833E-02  1.36667E-03 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  1.73627E+02  3.12783E+03 ];
CPU_USAGE                 (idx, 1)        = 0.97700 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.99298E-01 0.00010 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.68253E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  2.40000E-05 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 82 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  3.74315E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.31189E+01 ;
TOT_SF_RATE               (idx, 1)        =  4.06621E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  4.69276E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  3.41281E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  3.27386E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.27776E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  5.93980E+03 ;
INGESTION_TOXICITY        (idx, 1)        =  8.64601E+03 ;
ACTINIDE_INH_TOX          (idx, 1)        =  1.82506E+03 ;
ACTINIDE_ING_TOX          (idx, 1)        =  1.45197E+03 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  4.11474E+03 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  7.19404E+03 ;
SR90_ACTIVITY             (idx, 1)        =  7.10737E+07 ;
TE132_ACTIVITY            (idx, 1)        =  1.32878E+11 ;
I131_ACTIVITY             (idx, 1)        =  3.95563E+10 ;
I132_ACTIVITY             (idx, 1)        =  1.29663E+11 ;
CS134_ACTIVITY            (idx, 1)        =  2.24798E+04 ;
CS137_ACTIVITY            (idx, 1)        =  7.57238E+07 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  4.03655E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.97476E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.35714E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  5.51615E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.84565E+09 0.00099  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 3 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  1.00000E-01  1.00037E-01 ];
BURN_DAYS                 (idx, [1:  2])  = [  2.97885E+00  1.78731E+00 ];
FIMA                      (idx, [1:  3])  = [  1.05255E-04  1.67278E+18  1.58909E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.94446E-01 0.00234 ];
U235_FISS                 (idx, [1:   4]) = [  6.12887E+12 0.00117  9.44009E-01 0.00034 ];
U238_FISS                 (idx, [1:   4]) = [  3.57516E+11 0.00597  5.50438E-02 0.00572 ];
PU239_FISS                (idx, [1:   4]) = [  5.27124E+09 0.04600  8.12017E-04 0.04599 ];
U235_CAPT                 (idx, [1:   4]) = [  1.46211E+12 0.00281  2.80340E-01 0.00244 ];
U238_CAPT                 (idx, [1:   4]) = [  2.96137E+12 0.00224  5.67608E-01 0.00132 ];
PU239_CAPT                (idx, [1:   4]) = [  3.28562E+09 0.05731  6.30328E-04 0.05721 ];
PU240_CAPT                (idx, [1:   4]) = [  2.31165E+07 0.70645  4.58277E-06 0.70646 ];
XE135_CAPT                (idx, [1:   4]) = [  2.59499E+11 0.00659  4.97719E-02 0.00655 ];
SM149_CAPT                (idx, [1:   4]) = [  5.67537E+09 0.04515  1.09064E-03 0.04525 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000384 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51391E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000384 1.00151E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 445610 4.46129E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 554774 5.55385E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000384 1.00151E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 2.32831E-10 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.59783E+13 2.4E-05  1.59783E+13 2.4E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49692E+12 1.9E-06  6.49692E+12 1.9E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  5.22011E+12 0.00094  4.76097E+12 0.00103  4.59143E+11 0.00116 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.17170E+13 0.00042  1.12579E+13 0.00043  4.59143E+11 0.00116 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.16913E+13 0.00099  1.16913E+13 0.00099  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.91993E+14 0.00082  1.63386E+14 0.00094  3.28607E+14 0.00083 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.17170E+13 0.00042 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.14224E+14 0.00068 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27772E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27772E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.82463E+00 0.00078 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.53379E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.77474E-01 0.00087 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36030E+00 0.00082 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.36603E+00 0.00087 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.36603E+00 0.00087 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.45936E+00 2.5E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02479E+02 1.9E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.36615E+00 0.00094  1.35631E+00 0.00088  9.72266E-03 0.01612 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.36589E+00 0.00042 ];
COL_KEFF                  (idx, [1:   2]) = [  1.36735E+00 0.00099 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.36589E+00 0.00042 ];
ABS_KINF                  (idx, [1:   2]) = [  1.36589E+00 0.00042 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.69763E+01 0.00039 ];
IMP_ALF                   (idx, [1:   2]) = [  1.69802E+01 0.00018 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.57232E-07 0.00666 ];
IMP_EALF                  (idx, [1:   2]) = [  8.46440E-07 0.00304 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.97035E-01 0.00574 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.94914E-01 0.00245 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.29371E-03 0.01201  1.70812E-04 0.06449  7.94538E-04 0.03042  4.44413E-04 0.03923  9.86502E-04 0.02755  1.71240E-03 0.02123  5.64015E-04 0.03657  4.41926E-04 0.03856  1.79105E-04 0.06414 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.54656E-01 0.01829  4.68748E-03 0.05767  2.55757E-02 0.01459  3.07877E-02 0.02764  1.25326E-01 0.01111  2.90127E-01 0.00402  5.19860E-01 0.02377  1.17377E+00 0.02806  1.36497E+00 0.05670 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.08938E-03 0.01628  2.44039E-04 0.08970  1.08089E-03 0.04308  5.78533E-04 0.05443  1.26667E-03 0.04001  2.31495E-03 0.02814  7.61154E-04 0.05169  6.09419E-04 0.05390  2.33728E-04 0.09721 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.54565E-01 0.02544  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.8E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.0E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.21631E-05 0.00216  1.21582E-05 0.00217  1.29337E-05 0.02280 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.66087E-05 0.00190  1.66020E-05 0.00191  1.76604E-05 0.02273 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.12389E-03 0.01633  2.34772E-04 0.09042  1.04131E-03 0.04236  6.16357E-04 0.05440  1.29463E-03 0.03748  2.31708E-03 0.02881  7.54930E-04 0.05101  6.17309E-04 0.05382  2.47501E-04 0.08747 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.57982E-01 0.02664  1.24667E-02 0.0E+00  2.82917E-02 2.0E-09  4.25244E-02 8.1E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.5E-09  3.55460E+00 4.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.21325E-05 0.00448  1.21264E-05 0.00450  9.97808E-06 0.05160 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.65674E-05 0.00437  1.65590E-05 0.00439  1.36340E-05 0.05157 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.32791E-03 0.04338  1.27775E-04 0.26781  1.08361E-03 0.11101  8.38655E-04 0.13859  1.09894E-03 0.11700  2.46931E-03 0.08185  6.86794E-04 0.15234  6.26358E-04 0.14865  3.96473E-04 0.19461 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.10871E-01 0.06740  1.24667E-02 5.4E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 3.7E-09  2.92467E-01 6.1E-09  6.66488E-01 5.0E-09  1.63478E+00 0.0E+00  3.55460E+00 6.3E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.14709E-03 0.04169  1.39257E-04 0.26478  1.05722E-03 0.10745  8.06755E-04 0.13489  1.06710E-03 0.11453  2.43549E-03 0.07881  6.74024E-04 0.14912  6.19446E-04 0.14534  3.47798E-04 0.19458 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.97372E-01 0.06573  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 3.7E-09  2.92467E-01 6.1E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 6.8E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -6.11045E+02 0.04376 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.21949E-05 0.00123 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.66534E-05 0.00089 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.17723E-03 0.00842 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.89195E+02 0.00860 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.12458E-07 0.00120 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92894E-06 0.00103  2.92881E-06 0.00103  2.94827E-06 0.01268 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.92404E-05 0.00130  1.92420E-05 0.00131  1.90229E-05 0.01510 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.78162E-01 0.00087  5.76866E-01 0.00088  8.94739E-01 0.02149 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.27233E+01 0.02710 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.68428E+01 0.00057  2.93305E+01 0.00072 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.31882E+04 0.00643  5.41723E+04 0.00369  1.14170E+05 0.00237  1.25512E+05 0.00144  1.18418E+05 0.00184  1.32745E+05 0.00148  9.04990E+04 0.00227  8.11230E+04 0.00197  6.16534E+04 0.00226  5.02705E+04 0.00160  4.32095E+04 0.00199  3.90170E+04 0.00207  3.57708E+04 0.00222  3.41461E+04 0.00172  3.32418E+04 0.00173  2.85115E+04 0.00271  2.83130E+04 0.00225  2.78160E+04 0.00220  2.72207E+04 0.00228  5.25972E+04 0.00194  5.02221E+04 0.00117  3.57425E+04 0.00160  2.26870E+04 0.00236  2.56619E+04 0.00197  2.39311E+04 0.00249  2.18454E+04 0.00272  3.47104E+04 0.00148  8.02352E+03 0.00456  1.01295E+04 0.00266  9.27994E+03 0.00387  5.34382E+03 0.00489  9.33101E+03 0.00407  6.33069E+03 0.00374  5.29406E+03 0.00506  9.84179E+02 0.00907  9.78580E+02 0.01185  9.95704E+02 0.00889  1.03328E+03 0.00797  1.03738E+03 0.00806  1.02338E+03 0.00888  1.06153E+03 0.00819  9.83004E+02 0.01046  1.87979E+03 0.00735  3.01354E+03 0.00508  3.82909E+03 0.00602  1.00003E+04 0.00331  1.04452E+04 0.00370  1.09581E+04 0.00332  6.91330E+03 0.00437  4.79377E+03 0.00335  3.59761E+03 0.00629  4.05112E+03 0.00569  7.24806E+03 0.00314  8.99853E+03 0.00260  1.60722E+04 0.00196  2.17651E+04 0.00222  2.91008E+04 0.00226  1.71119E+04 0.00260  1.16506E+04 0.00314  8.13966E+03 0.00309  7.14239E+03 0.00282  6.83404E+03 0.00376  5.61518E+03 0.00411  3.73452E+03 0.00525  3.39466E+03 0.00404  2.96992E+03 0.00386  2.46099E+03 0.00527  1.89893E+03 0.00573  1.23486E+03 0.00595  4.25893E+02 0.00987 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.36689E+00 0.00118 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.34607E+14 0.00109  5.74483E+13 0.00076 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.46121E-01 0.00022  1.33924E+00 0.00039 ];
INF_CAPT                  (idx, [1:   4]) = [  7.56591E-03 0.00101  3.36476E-02 0.00055 ];
INF_ABS                   (idx, [1:   4]) = [  1.14231E-02 0.00076  1.17612E-01 0.00068 ];
INF_FISS                  (idx, [1:   4]) = [  3.85718E-03 0.00089  8.39645E-02 0.00073 ];
INF_NSF                   (idx, [1:   4]) = [  9.73820E-03 0.00088  2.04592E-01 0.00073 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52469E+00 7.8E-05  2.43664E+00 2.0E-07 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03064E+02 6.8E-06  2.02276E+02 3.3E-08 ];
INF_INVV                  (idx, [1:   4]) = [  5.50266E-08 0.00110  2.26011E-06 0.00041 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34715E-01 0.00024  1.22164E+00 0.00047 ];
INF_SCATT1                (idx, [1:   4]) = [  2.36133E-01 0.00047  3.31208E-01 0.00105 ];
INF_SCATT2                (idx, [1:   4]) = [  9.36348E-02 0.00069  8.45704E-02 0.00265 ];
INF_SCATT3                (idx, [1:   4]) = [  7.42497E-03 0.00545  2.53053E-02 0.00871 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.47334E-03 0.00414 -5.59042E-03 0.03002 ];
INF_SCATT5                (idx, [1:   4]) = [  4.74494E-04 0.07830  4.30029E-03 0.03187 ];
INF_SCATT6                (idx, [1:   4]) = [  5.00716E-03 0.00661 -1.23207E-02 0.01175 ];
INF_SCATT7                (idx, [1:   4]) = [  7.37885E-04 0.04795 -4.53620E-04 0.23096 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34756E-01 0.00024  1.22164E+00 0.00047 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.36135E-01 0.00047  3.31208E-01 0.00105 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.36336E-02 0.00069  8.45704E-02 0.00265 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.42554E-03 0.00546  2.53053E-02 0.00871 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.47388E-03 0.00414 -5.59042E-03 0.03002 ];
INF_SCATTP5               (idx, [1:   4]) = [  4.74248E-04 0.07834  4.30029E-03 0.03187 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.00739E-03 0.00661 -1.23207E-02 0.01175 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.37270E-04 0.04806 -4.53620E-04 0.23096 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30740E-01 0.00035  8.85386E-01 0.00044 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44463E+00 0.00035  3.76485E-01 0.00044 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13824E-02 0.00079  1.17612E-01 0.00068 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72439E-02 0.00045  1.19893E-01 0.00100 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18877E-01 0.00023  1.58380E-02 0.00082  2.29566E-03 0.01047  1.21935E+00 0.00047 ];
INF_S1                    (idx, [1:   8]) = [  2.31550E-01 0.00047  4.58372E-03 0.00169  9.03356E-04 0.01924  3.30304E-01 0.00107 ];
INF_S2                    (idx, [1:   8]) = [  9.50159E-02 0.00069 -1.38105E-03 0.00427  4.96102E-04 0.02559  8.40743E-02 0.00271 ];
INF_S3                    (idx, [1:   8]) = [  9.04922E-03 0.00441 -1.62424E-03 0.00325  1.63360E-04 0.04832  2.51420E-02 0.00870 ];
INF_S4                    (idx, [1:   8]) = [ -8.94664E-03 0.00422 -5.26709E-04 0.01061 -1.07641E-05 0.60052 -5.57966E-03 0.02980 ];
INF_S5                    (idx, [1:   8]) = [  4.52127E-04 0.08155  2.23666E-05 0.19597 -7.24567E-05 0.08503  4.37275E-03 0.03162 ];
INF_S6                    (idx, [1:   8]) = [  5.13682E-03 0.00642 -1.29665E-04 0.04492 -1.02944E-04 0.04568 -1.22177E-02 0.01192 ];
INF_S7                    (idx, [1:   8]) = [  8.93385E-04 0.03956 -1.55499E-04 0.03792 -9.12864E-05 0.06230 -3.62334E-04 0.29516 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18918E-01 0.00023  1.58380E-02 0.00082  2.29566E-03 0.01047  1.21935E+00 0.00047 ];
INF_SP1                   (idx, [1:   8]) = [  2.31551E-01 0.00046  4.58372E-03 0.00169  9.03356E-04 0.01924  3.30304E-01 0.00107 ];
INF_SP2                   (idx, [1:   8]) = [  9.50147E-02 0.00068 -1.38105E-03 0.00427  4.96102E-04 0.02559  8.40743E-02 0.00271 ];
INF_SP3                   (idx, [1:   8]) = [  9.04978E-03 0.00442 -1.62424E-03 0.00325  1.63360E-04 0.04832  2.51420E-02 0.00870 ];
INF_SP4                   (idx, [1:   8]) = [ -8.94717E-03 0.00422 -5.26709E-04 0.01061 -1.07641E-05 0.60052 -5.57966E-03 0.02980 ];
INF_SP5                   (idx, [1:   8]) = [  4.51881E-04 0.08159  2.23666E-05 0.19597 -7.24567E-05 0.08503  4.37275E-03 0.03162 ];
INF_SP6                   (idx, [1:   8]) = [  5.13705E-03 0.00643 -1.29665E-04 0.04492 -1.02944E-04 0.04568 -1.22177E-02 0.01192 ];
INF_SP7                   (idx, [1:   8]) = [  8.92769E-04 0.03962 -1.55499E-04 0.03792 -9.12864E-05 0.06230 -3.62334E-04 0.29516 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42742E-01 0.00143  7.91558E-01 0.00617 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.42591E-01 0.00207  7.92326E-01 0.00807 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.43236E-01 0.00220  7.86875E-01 0.01139 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.42435E-01 0.00175  7.99824E-01 0.01181 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37327E+00 0.00142  4.21495E-01 0.00616 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37420E+00 0.00208  4.21364E-01 0.00812 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37057E+00 0.00218  4.24953E-01 0.01152 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37504E+00 0.00176  4.18167E-01 0.01191 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.08938E-03 0.01628  2.44039E-04 0.08970  1.08089E-03 0.04308  5.78533E-04 0.05443  1.26667E-03 0.04001  2.31495E-03 0.02814  7.61154E-04 0.05169  6.09419E-04 0.05390  2.33728E-04 0.09721 ];
LAMBDA                    (idx, [1:  18]) = [  4.54565E-01 0.02544  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.8E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.0E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Wed Feb  7 22:42:48 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 0.0E+00  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.68331E-02 0.00142  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73167E-01 3.9E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.34367E-01 0.00014  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98773E-01 7.0E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98272E-01 9.7E-06 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.34061E-01 0.00014  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.82433E+00 0.00050  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.68567E+01 0.00053  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.68567E+01 0.00053  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.68516E+00 0.00072  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.20747E-01 0.00154  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000540 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00108E+03 0.00131 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00108E+03 0.00131 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.29314E+02 ;
RUNNING_TIME              (idx, 1)        =  2.33405E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.32127E+00  1.88650E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  2.30122E+02  2.97644E+01  2.96350E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  9.19667E-02  5.76667E-03 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  3.67167E-02  1.35000E-03 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.33405E+02  3.25397E+03 ];
CPU_USAGE                 (idx, 1)        = 0.98247 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.99096E-01 9.8E-05 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.76105E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  2.40000E-05 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 88 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  3.91207E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.32850E+01 ;
TOT_SF_RATE               (idx, 1)        =  4.10942E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  5.43464E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  3.94256E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  3.36860E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.28907E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  8.46502E+03 ;
INGESTION_TOXICITY        (idx, 1)        =  1.10578E+04 ;
ACTINIDE_INH_TOX          (idx, 1)        =  2.62288E+03 ;
ACTINIDE_ING_TOX          (idx, 1)        =  2.03948E+03 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  5.84214E+03 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  9.01832E+03 ;
SR90_ACTIVITY             (idx, 1)        =  1.42127E+08 ;
TE132_ACTIVITY            (idx, 1)        =  2.02735E+11 ;
I131_ACTIVITY             (idx, 1)        =  7.30186E+10 ;
I132_ACTIVITY             (idx, 1)        =  2.01685E+11 ;
CS134_ACTIVITY            (idx, 1)        =  2.04551E+05 ;
CS137_ACTIVITY            (idx, 1)        =  1.51534E+08 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  4.27251E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.96708E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.40636E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  6.11741E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.86633E+09 0.00099  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 4 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  2.00000E-01  2.00075E-01 ];
BURN_DAYS                 (idx, [1:  2])  = [  5.95770E+00  2.97885E+00 ];
FIMA                      (idx, [1:  3])  = [  2.10502E-04  3.34542E+18  1.58892E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.94433E-01 0.00236 ];
U235_FISS                 (idx, [1:   4]) = [  6.12500E+12 0.00103  9.41865E-01 0.00034 ];
U238_FISS                 (idx, [1:   4]) = [  3.59140E+11 0.00593  5.51904E-02 0.00562 ];
PU239_FISS                (idx, [1:   4]) = [  1.79747E+10 0.02636  2.76543E-03 0.02636 ];
PU240_FISS                (idx, [1:   4]) = [  1.10843E+07 1.00000  1.68776E-06 1.00000 ];
U235_CAPT                 (idx, [1:   4]) = [  1.45997E+12 0.00272  2.78331E-01 0.00242 ];
U238_CAPT                 (idx, [1:   4]) = [  2.96641E+12 0.00230  5.65249E-01 0.00135 ];
PU239_CAPT                (idx, [1:   4]) = [  1.00303E+10 0.03398  1.91038E-03 0.03395 ];
PU240_CAPT                (idx, [1:   4]) = [  1.16856E+08 0.31343  2.14954E-05 0.31372 ];
XE135_CAPT                (idx, [1:   4]) = [  2.62896E+11 0.00647  5.01236E-02 0.00639 ];
SM149_CAPT                (idx, [1:   4]) = [  1.78909E+10 0.02542  3.41429E-03 0.02547 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000540 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.51674E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000540 1.00152E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 446694 4.47143E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 553846 5.54374E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000540 1.00152E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -2.32831E-09 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.59830E+13 2.2E-05  1.59830E+13 2.2E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49656E+12 2.0E-06  6.49656E+12 2.0E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  5.24635E+12 0.00093  4.78637E+12 0.00101  4.59981E+11 0.00119 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.17429E+13 0.00042  1.12829E+13 0.00043  4.59981E+11 0.00119 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.17327E+13 0.00099  1.17327E+13 0.00099  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.93634E+14 0.00087  1.63860E+14 0.00096  3.29774E+14 0.00089 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.17429E+13 0.00042 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.15123E+14 0.00071 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27706E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27706E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.81903E+00 0.00076 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.53718E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.77505E-01 0.00084 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36190E+00 0.00080 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.36403E+00 0.00086 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.36403E+00 0.00086 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46023E+00 2.4E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02490E+02 2.0E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.36386E+00 0.00092  1.35443E+00 0.00086  9.59926E-03 0.01541 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.36326E+00 0.00041 ];
COL_KEFF                  (idx, [1:   2]) = [  1.36292E+00 0.00098 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.36326E+00 0.00041 ];
ABS_KINF                  (idx, [1:   2]) = [  1.36326E+00 0.00041 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.69648E+01 0.00039 ];
IMP_ALF                   (idx, [1:   2]) = [  1.69711E+01 0.00019 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.66946E-07 0.00665 ];
IMP_EALF                  (idx, [1:   2]) = [  8.54346E-07 0.00315 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.97526E-01 0.00578 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.95564E-01 0.00231 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.35352E-03 0.01135  1.59344E-04 0.06885  8.04245E-04 0.02978  4.70752E-04 0.04035  1.01861E-03 0.02671  1.66964E-03 0.02057  5.78824E-04 0.03670  4.77347E-04 0.03907  1.74756E-04 0.06150 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.61828E-01 0.01885  4.33841E-03 0.06128  2.50099E-02 0.01622  3.06176E-02 0.02792  1.25326E-01 0.01111  2.88373E-01 0.00533  5.26525E-01 0.02308  1.20974E+00 0.02654  1.42895E+00 0.05460 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.18567E-03 0.01575  2.14351E-04 0.09745  1.07116E-03 0.04043  6.22042E-04 0.05426  1.35564E-03 0.03577  2.20148E-03 0.03033  8.28427E-04 0.05328  6.79966E-04 0.05339  2.12612E-04 0.08547 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.59331E-01 0.02361  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.21397E-05 0.00207  1.21325E-05 0.00208  1.30652E-05 0.02169 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.65492E-05 0.00181  1.65393E-05 0.00182  1.78051E-05 0.02159 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.07318E-03 0.01619  1.93531E-04 0.09499  1.06742E-03 0.04221  6.61239E-04 0.05209  1.40905E-03 0.03774  2.11679E-03 0.03201  7.64480E-04 0.04985  6.19292E-04 0.05434  2.41383E-04 0.08588 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.67114E-01 0.02708  1.24667E-02 0.0E+00  2.82917E-02 1.7E-09  4.25244E-02 8.0E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.6E-09  3.55460E+00 4.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.21824E-05 0.00437  1.21775E-05 0.00438  1.02610E-05 0.04915 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.66073E-05 0.00425  1.66007E-05 0.00426  1.39681E-05 0.04898 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.20870E-03 0.04703  1.97790E-04 0.29690  1.04553E-03 0.12752  8.23543E-04 0.15122  1.47329E-03 0.10607  2.17009E-03 0.08650  6.02336E-04 0.15400  6.14985E-04 0.16490  2.81136E-04 0.23724 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.97923E-01 0.07283  1.24667E-02 4.7E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.7E-09  2.92467E-01 5.8E-09  6.66488E-01 5.0E-09  1.63478E+00 0.0E+00  3.55460E+00 0.0E+00 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.24216E-03 0.04735  2.00045E-04 0.29318  1.04896E-03 0.12522  8.17566E-04 0.15253  1.43749E-03 0.10475  2.19211E-03 0.08427  6.39134E-04 0.14939  6.31436E-04 0.16348  2.75424E-04 0.23788 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.94065E-01 0.07150  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.9E-09  2.92467E-01 5.9E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 6.6E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -5.97940E+02 0.04730 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.21760E-05 0.00122 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.65995E-05 0.00082 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.00690E-03 0.00852 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.75863E+02 0.00863 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.11931E-07 0.00115 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92384E-06 0.00103  2.92400E-06 0.00103  2.91599E-06 0.01262 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.91745E-05 0.00124  1.91769E-05 0.00124  1.89739E-05 0.01557 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.78161E-01 0.00084  5.76823E-01 0.00084  8.87720E-01 0.02028 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.24134E+01 0.02547 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.68237E+01 0.00053  2.92784E+01 0.00068 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.30176E+04 0.00623  5.43226E+04 0.00284  1.14453E+05 0.00224  1.25704E+05 0.00198  1.18101E+05 0.00149  1.32634E+05 0.00146  9.02479E+04 0.00159  8.11948E+04 0.00197  6.17079E+04 0.00184  5.03862E+04 0.00203  4.32783E+04 0.00206  3.90877E+04 0.00215  3.59068E+04 0.00189  3.41302E+04 0.00174  3.30754E+04 0.00182  2.85539E+04 0.00158  2.81382E+04 0.00279  2.78038E+04 0.00243  2.72511E+04 0.00180  5.26282E+04 0.00167  5.01024E+04 0.00184  3.58147E+04 0.00182  2.26453E+04 0.00237  2.57659E+04 0.00214  2.40142E+04 0.00210  2.17220E+04 0.00247  3.47766E+04 0.00204  8.08347E+03 0.00319  1.01595E+04 0.00399  9.24977E+03 0.00358  5.33468E+03 0.00412  9.29870E+03 0.00249  6.26295E+03 0.00618  5.24919E+03 0.00568  9.76572E+02 0.00941  9.80612E+02 0.00890  1.01810E+03 0.01054  1.04246E+03 0.00973  1.02252E+03 0.00985  1.00268E+03 0.00764  1.06353E+03 0.00852  9.92637E+02 0.00966  1.88231E+03 0.00874  3.01497E+03 0.00427  3.82476E+03 0.00523  1.00386E+04 0.00387  1.04789E+04 0.00339  1.09215E+04 0.00380  6.96539E+03 0.00378  4.80455E+03 0.00474  3.62335E+03 0.00528  4.03336E+03 0.00396  7.22315E+03 0.00281  8.95839E+03 0.00424  1.60759E+04 0.00230  2.17128E+04 0.00216  2.89480E+04 0.00182  1.70886E+04 0.00203  1.16265E+04 0.00267  8.14714E+03 0.00362  7.12960E+03 0.00293  6.76715E+03 0.00262  5.58495E+03 0.00270  3.69886E+03 0.00400  3.36160E+03 0.00480  2.96677E+03 0.00350  2.45484E+03 0.00430  1.91515E+03 0.00436  1.24171E+03 0.00571  4.31794E+02 0.00978 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.36322E+00 0.00087 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.36150E+14 0.00092  5.75360E+13 0.00058 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.46031E-01 0.00024  1.33880E+00 0.00040 ];
INF_CAPT                  (idx, [1:   4]) = [  7.55828E-03 0.00116  3.39037E-02 0.00046 ];
INF_ABS                   (idx, [1:   4]) = [  1.14096E-02 0.00080  1.17676E-01 0.00051 ];
INF_FISS                  (idx, [1:   4]) = [  3.85128E-03 0.00076  8.37720E-02 0.00053 ];
INF_NSF                   (idx, [1:   4]) = [  9.72473E-03 0.00075  2.04204E-01 0.00053 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52506E+00 7.7E-05  2.43762E+00 5.6E-07 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03069E+02 8.0E-06  2.02289E+02 9.1E-08 ];
INF_INVV                  (idx, [1:   4]) = [  5.50068E-08 0.00080  2.25965E-06 0.00036 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34628E-01 0.00026  1.22106E+00 0.00045 ];
INF_SCATT1                (idx, [1:   4]) = [  2.36111E-01 0.00047  3.30985E-01 0.00125 ];
INF_SCATT2                (idx, [1:   4]) = [  9.35862E-02 0.00090  8.47179E-02 0.00268 ];
INF_SCATT3                (idx, [1:   4]) = [  7.39454E-03 0.00686  2.56804E-02 0.00833 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.40175E-03 0.00358 -5.30204E-03 0.02603 ];
INF_SCATT5                (idx, [1:   4]) = [  4.66293E-04 0.06401  4.42398E-03 0.02164 ];
INF_SCATT6                (idx, [1:   4]) = [  5.07502E-03 0.00760 -1.22909E-02 0.01143 ];
INF_SCATT7                (idx, [1:   4]) = [  7.63761E-04 0.04787 -4.92292E-04 0.28476 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34669E-01 0.00026  1.22106E+00 0.00045 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.36111E-01 0.00047  3.30985E-01 0.00125 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.35863E-02 0.00090  8.47179E-02 0.00268 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.39509E-03 0.00688  2.56804E-02 0.00833 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.40189E-03 0.00359 -5.30204E-03 0.02603 ];
INF_SCATTP5               (idx, [1:   4]) = [  4.65996E-04 0.06420  4.42398E-03 0.02164 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.07457E-03 0.00761 -1.22909E-02 0.01143 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.63776E-04 0.04774 -4.92292E-04 0.28476 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30773E-01 0.00058  8.85750E-01 0.00048 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44443E+00 0.00058  3.76331E-01 0.00048 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13688E-02 0.00077  1.17676E-01 0.00051 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72438E-02 0.00052  1.20055E-01 0.00093 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18788E-01 0.00026  1.58404E-02 0.00100  2.31339E-03 0.00948  1.21875E+00 0.00046 ];
INF_S1                    (idx, [1:   8]) = [  2.31532E-01 0.00047  4.57951E-03 0.00193  9.23396E-04 0.01692  3.30061E-01 0.00124 ];
INF_S2                    (idx, [1:   8]) = [  9.49785E-02 0.00088 -1.39224E-03 0.00550  4.90493E-04 0.02402  8.42274E-02 0.00268 ];
INF_S3                    (idx, [1:   8]) = [  9.01362E-03 0.00538 -1.61909E-03 0.00541  1.64067E-04 0.05376  2.55163E-02 0.00837 ];
INF_S4                    (idx, [1:   8]) = [ -8.87506E-03 0.00402 -5.26694E-04 0.01416 -5.30123E-07 1.00000 -5.30151E-03 0.02593 ];
INF_S5                    (idx, [1:   8]) = [  4.49093E-04 0.06676  1.71997E-05 0.33170 -7.73835E-05 0.06377  4.50136E-03 0.02091 ];
INF_S6                    (idx, [1:   8]) = [  5.19533E-03 0.00721 -1.20304E-04 0.04900 -9.90210E-05 0.06717 -1.21918E-02 0.01165 ];
INF_S7                    (idx, [1:   8]) = [  9.15348E-04 0.03958 -1.51586E-04 0.03037 -9.25387E-05 0.06860 -3.99754E-04 0.35126 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18828E-01 0.00026  1.58404E-02 0.00100  2.31339E-03 0.00948  1.21875E+00 0.00046 ];
INF_SP1                   (idx, [1:   8]) = [  2.31532E-01 0.00047  4.57951E-03 0.00193  9.23396E-04 0.01692  3.30061E-01 0.00124 ];
INF_SP2                   (idx, [1:   8]) = [  9.49786E-02 0.00088 -1.39224E-03 0.00550  4.90493E-04 0.02402  8.42274E-02 0.00268 ];
INF_SP3                   (idx, [1:   8]) = [  9.01418E-03 0.00540 -1.61909E-03 0.00541  1.64067E-04 0.05376  2.55163E-02 0.00837 ];
INF_SP4                   (idx, [1:   8]) = [ -8.87520E-03 0.00403 -5.26694E-04 0.01416 -5.30123E-07 1.00000 -5.30151E-03 0.02593 ];
INF_SP5                   (idx, [1:   8]) = [  4.48797E-04 0.06699  1.71997E-05 0.33170 -7.73835E-05 0.06377  4.50136E-03 0.02091 ];
INF_SP6                   (idx, [1:   8]) = [  5.19487E-03 0.00721 -1.20304E-04 0.04900 -9.90210E-05 0.06717 -1.21918E-02 0.01165 ];
INF_SP7                   (idx, [1:   8]) = [  9.15363E-04 0.03946 -1.51586E-04 0.03037 -9.25387E-05 0.06860 -3.99754E-04 0.35126 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42830E-01 0.00119  7.94896E-01 0.00350 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.43009E-01 0.00160  7.95644E-01 0.00761 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.43340E-01 0.00206  8.03014E-01 0.00863 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.42173E-01 0.00163  7.89413E-01 0.00881 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37275E+00 0.00119  4.19466E-01 0.00352 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37177E+00 0.00161  4.19521E-01 0.00749 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.36997E+00 0.00207  4.15836E-01 0.00853 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37651E+00 0.00162  4.23042E-01 0.00881 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.18567E-03 0.01575  2.14351E-04 0.09745  1.07116E-03 0.04043  6.22042E-04 0.05426  1.35564E-03 0.03577  2.20148E-03 0.03033  8.28427E-04 0.05328  6.79966E-04 0.05339  2.12612E-04 0.08547 ];
LAMBDA                    (idx, [1:  18]) = [  4.59331E-01 0.02361  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Wed Feb  7 23:43:53 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.6E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.68360E-02 0.00141  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73164E-01 3.9E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.34479E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98773E-01 7.1E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98271E-01 1.0E-05 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.34171E-01 0.00013  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.82399E+00 0.00049  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.68381E+01 0.00054  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.68381E+01 0.00054  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.67291E+00 0.00071  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.20398E-01 0.00155  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000304 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00061E+03 0.00137 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00061E+03 0.00137 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  2.90016E+02 ;
RUNNING_TIME              (idx, 1)        =  2.94493E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  1.72117E+00  2.00633E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  2.90778E+02  3.10407E+01  2.96150E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  1.21933E-01  2.35833E-02 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  4.73333E-02  1.58333E-03 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  2.94493E+02  3.28589E+03 ];
CPU_USAGE                 (idx, 1)        = 0.98480 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.99038E-01 0.00011 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.79942E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.10000E-06 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 90 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  4.02509E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.33999E+01 ;
TOT_SF_RATE               (idx, 1)        =  4.25224E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  5.82137E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  4.21816E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  3.44295E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.29780E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  1.09052E+04 ;
INGESTION_TOXICITY        (idx, 1)        =  1.30443E+04 ;
ACTINIDE_INH_TOX          (idx, 1)        =  3.12611E+03 ;
ACTINIDE_ING_TOX          (idx, 1)        =  2.34774E+03 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  7.77905E+03 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  1.06965E+04 ;
SR90_ACTIVITY             (idx, 1)        =  2.48470E+08 ;
TE132_ACTIVITY            (idx, 1)        =  2.50662E+11 ;
I131_ACTIVITY             (idx, 1)        =  1.10250E+11 ;
I132_ACTIVITY             (idx, 1)        =  2.51116E+11 ;
CS134_ACTIVITY            (idx, 1)        =  1.17817E+06 ;
CS137_ACTIVITY            (idx, 1)        =  2.65231E+08 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  4.41928E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.95299E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.49929E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  6.47591E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.88433E+09 0.00099  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 5 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  3.50000E-01  3.50132E-01 ];
BURN_DAYS                 (idx, [1:  2])  = [  1.04260E+01  4.46828E+00 ];
FIMA                      (idx, [1:  3])  = [  3.68362E-04  5.85422E+18  1.58867E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.95569E-01 0.00231 ];
U235_FISS                 (idx, [1:   4]) = [  6.09257E+12 0.00112  9.38205E-01 0.00032 ];
U238_FISS                 (idx, [1:   4]) = [  3.60609E+11 0.00555  5.54965E-02 0.00522 ];
PU239_FISS                (idx, [1:   4]) = [  3.98415E+10 0.01665  6.13283E-03 0.01657 ];
PU240_FISS                (idx, [1:   4]) = [  1.23678E+07 1.00000  1.90658E-06 1.00000 ];
PU241_FISS                (idx, [1:   4]) = [  1.21124E+07 1.00000  1.82482E-06 1.00000 ];
U235_CAPT                 (idx, [1:   4]) = [  1.46497E+12 0.00271  2.76849E-01 0.00238 ];
U238_CAPT                 (idx, [1:   4]) = [  2.97722E+12 0.00227  5.62391E-01 0.00131 ];
PU239_CAPT                (idx, [1:   4]) = [  2.36543E+10 0.02209  4.46404E-03 0.02196 ];
PU240_CAPT                (idx, [1:   4]) = [  4.60312E+08 0.16181  8.65307E-05 0.16250 ];
XE135_CAPT                (idx, [1:   4]) = [  2.60492E+11 0.00668  4.92101E-02 0.00644 ];
SM149_CAPT                (idx, [1:   4]) = [  3.29831E+10 0.01919  6.23865E-03 0.01924 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000304 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.57362E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000304 1.00157E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 449076 4.49692E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 551228 5.51882E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000304 1.00157E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 4.65661E-10 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.59924E+13 2.3E-05  1.59924E+13 2.3E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49588E+12 2.0E-06  6.49588E+12 2.0E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  5.28835E+12 0.00093  4.82874E+12 0.00100  4.59616E+11 0.00123 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.17842E+13 0.00042  1.13246E+13 0.00042  4.59616E+11 0.00123 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.17687E+13 0.00099  1.17687E+13 0.00099  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.95052E+14 0.00084  1.64514E+14 0.00091  3.30537E+14 0.00087 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.17842E+13 0.00042 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.15883E+14 0.00072 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27609E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27609E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.81067E+00 0.00080 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.54098E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.77200E-01 0.00087 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36304E+00 0.00078 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.35872E+00 0.00091 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.35872E+00 0.00091 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46193E+00 2.4E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02512E+02 2.0E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.35843E+00 0.00098  1.34885E+00 0.00092  9.86987E-03 0.01589 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.35931E+00 0.00042 ];
COL_KEFF                  (idx, [1:   2]) = [  1.35956E+00 0.00099 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.35931E+00 0.00042 ];
ABS_KINF                  (idx, [1:   2]) = [  1.35931E+00 0.00042 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.69540E+01 0.00038 ];
IMP_ALF                   (idx, [1:   2]) = [  1.69617E+01 0.00017 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.75879E-07 0.00644 ];
IMP_EALF                  (idx, [1:   2]) = [  8.62164E-07 0.00294 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.98345E-01 0.00558 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.97358E-01 0.00233 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.41181E-03 0.01159  1.59225E-04 0.07077  7.76020E-04 0.03123  4.91284E-04 0.03861  1.01695E-03 0.02761  1.73522E-03 0.02069  5.40645E-04 0.03658  4.83211E-04 0.03936  2.09251E-04 0.05905 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.81255E-01 0.01932  4.21374E-03 0.06265  2.48401E-02 0.01669  3.17232E-02 0.02612  1.24261E-01 0.01190  2.91297E-01 0.00284  5.09197E-01 0.02488  1.18031E+00 0.02778  1.55691E+00 0.05071 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.08126E-03 0.01617  2.16241E-04 0.09286  9.67851E-04 0.04562  6.48742E-04 0.05413  1.37145E-03 0.04156  2.30800E-03 0.02826  7.06620E-04 0.05239  6.05925E-04 0.05564  2.56431E-04 0.08624 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.71813E-01 0.02617  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.6E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.21250E-05 0.00211  1.21184E-05 0.00212  1.30063E-05 0.02336 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.64627E-05 0.00185  1.64538E-05 0.00186  1.76515E-05 0.02323 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.29109E-03 0.01615  1.91486E-04 0.10055  1.03089E-03 0.04320  6.82456E-04 0.05177  1.36119E-03 0.03826  2.38443E-03 0.02843  7.26624E-04 0.05235  6.74587E-04 0.05141  2.39425E-04 0.09050 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.71948E-01 0.02890  1.24667E-02 0.0E+00  2.82917E-02 1.9E-09  4.25244E-02 7.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.5E-09  3.55460E+00 4.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.21144E-05 0.00430  1.21096E-05 0.00433  1.04028E-05 0.04868 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.64512E-05 0.00425  1.64446E-05 0.00428  1.41438E-05 0.04865 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.41250E-03 0.04354  2.51630E-04 0.29152  9.88898E-04 0.12059  7.11868E-04 0.13610  1.28006E-03 0.11323  2.55650E-03 0.07683  7.68069E-04 0.15558  6.38677E-04 0.13419  2.16791E-04 0.24049 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.78821E-01 0.06731  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.4E-09  2.92467E-01 6.2E-09  6.66488E-01 5.5E-09  1.63478E+00 0.0E+00  3.55460E+00 5.4E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.35448E-03 0.04251  2.36266E-04 0.27521  9.68299E-04 0.11797  7.24588E-04 0.13316  1.27174E-03 0.10995  2.52346E-03 0.07361  7.66586E-04 0.14891  6.46637E-04 0.13380  2.16896E-04 0.22917 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.75336E-01 0.06568  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.2E-09  2.92467E-01 6.2E-09  6.66488E-01 5.0E-09  1.63478E+00 0.0E+00  3.55460E+00 4.6E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -6.20922E+02 0.04487 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.21357E-05 0.00130 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.64774E-05 0.00083 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.31836E-03 0.00840 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -6.03669E+02 0.00853 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.10698E-07 0.00117 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92695E-06 0.00100  2.92689E-06 0.00100  2.92955E-06 0.01209 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.91181E-05 0.00132  1.91196E-05 0.00132  1.88364E-05 0.01518 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.77851E-01 0.00087  5.76579E-01 0.00087  8.82784E-01 0.02209 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.23797E+01 0.02705 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.68051E+01 0.00054  2.92534E+01 0.00068 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.33257E+04 0.00500  5.43416E+04 0.00325  1.14518E+05 0.00221  1.25909E+05 0.00189  1.18116E+05 0.00149  1.32870E+05 0.00152  9.04224E+04 0.00147  8.10742E+04 0.00169  6.18012E+04 0.00202  5.03143E+04 0.00135  4.32167E+04 0.00204  3.90662E+04 0.00159  3.59345E+04 0.00174  3.40067E+04 0.00178  3.30663E+04 0.00227  2.85989E+04 0.00246  2.81998E+04 0.00175  2.77859E+04 0.00245  2.72650E+04 0.00181  5.27510E+04 0.00132  5.01471E+04 0.00177  3.57462E+04 0.00213  2.27658E+04 0.00231  2.56158E+04 0.00181  2.39124E+04 0.00240  2.17077E+04 0.00310  3.47988E+04 0.00198  8.03063E+03 0.00467  1.00914E+04 0.00430  9.19755E+03 0.00295  5.33127E+03 0.00437  9.37676E+03 0.00348  6.27340E+03 0.00423  5.27944E+03 0.00572  1.00746E+03 0.00914  9.65348E+02 0.00756  9.94293E+02 0.00833  1.04869E+03 0.01080  1.03592E+03 0.00837  1.00097E+03 0.00959  1.05372E+03 0.00886  9.89432E+02 0.00941  1.90598E+03 0.00765  3.00044E+03 0.00554  3.79011E+03 0.00576  9.97927E+03 0.00320  1.03592E+04 0.00306  1.08999E+04 0.00377  6.85984E+03 0.00417  4.74866E+03 0.00399  3.57122E+03 0.00399  4.07331E+03 0.00564  7.18404E+03 0.00285  8.94124E+03 0.00315  1.59846E+04 0.00342  2.16363E+04 0.00227  2.88705E+04 0.00265  1.69872E+04 0.00230  1.15101E+04 0.00275  8.07259E+03 0.00245  7.11533E+03 0.00352  6.79525E+03 0.00385  5.58956E+03 0.00293  3.70532E+03 0.00328  3.36058E+03 0.00347  2.93547E+03 0.00438  2.44551E+03 0.00472  1.88447E+03 0.00490  1.23870E+03 0.00614  4.26432E+02 0.00825 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.35976E+00 0.00097 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.37685E+14 0.00110  5.74282E+13 0.00090 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.45862E-01 0.00027  1.34125E+00 0.00045 ];
INF_CAPT                  (idx, [1:   4]) = [  7.56892E-03 0.00107  3.44168E-02 0.00063 ];
INF_ABS                   (idx, [1:   4]) = [  1.14139E-02 0.00071  1.18280E-01 0.00078 ];
INF_FISS                  (idx, [1:   4]) = [  3.84495E-03 0.00094  8.38632E-02 0.00085 ];
INF_NSF                   (idx, [1:   4]) = [  9.71343E-03 0.00094  2.04578E-01 0.00085 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52628E+00 6.1E-05  2.43943E+00 1.6E-06 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03080E+02 5.6E-06  2.02313E+02 2.6E-07 ];
INF_INVV                  (idx, [1:   4]) = [  5.49190E-08 0.00095  2.26016E-06 0.00039 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34454E-01 0.00028  1.22298E+00 0.00055 ];
INF_SCATT1                (idx, [1:   4]) = [  2.36009E-01 0.00052  3.31749E-01 0.00124 ];
INF_SCATT2                (idx, [1:   4]) = [  9.35457E-02 0.00078  8.51180E-02 0.00296 ];
INF_SCATT3                (idx, [1:   4]) = [  7.44581E-03 0.00711  2.59802E-02 0.00784 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.39034E-03 0.00443 -5.23476E-03 0.02994 ];
INF_SCATT5                (idx, [1:   4]) = [  5.27780E-04 0.06879  4.32613E-03 0.03782 ];
INF_SCATT6                (idx, [1:   4]) = [  5.08252E-03 0.00805 -1.24625E-02 0.01256 ];
INF_SCATT7                (idx, [1:   4]) = [  7.10141E-04 0.05509 -7.30269E-04 0.20010 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34496E-01 0.00028  1.22298E+00 0.00055 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.36010E-01 0.00052  3.31749E-01 0.00124 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.35460E-02 0.00078  8.51180E-02 0.00296 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.44597E-03 0.00711  2.59802E-02 0.00784 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.39070E-03 0.00443 -5.23476E-03 0.02994 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.27285E-04 0.06891  4.32613E-03 0.03782 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.08254E-03 0.00805 -1.24625E-02 0.01256 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.10067E-04 0.05530 -7.30269E-04 0.20010 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30658E-01 0.00044  8.87081E-01 0.00052 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44515E+00 0.00044  3.75767E-01 0.00052 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13716E-02 0.00071  1.18280E-01 0.00078 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72292E-02 0.00052  1.20571E-01 0.00122 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18632E-01 0.00027  1.58212E-02 0.00088  2.29507E-03 0.00628  1.22068E+00 0.00056 ];
INF_S1                    (idx, [1:   8]) = [  2.31438E-01 0.00050  4.57101E-03 0.00225  9.22498E-04 0.01257  3.30826E-01 0.00124 ];
INF_S2                    (idx, [1:   8]) = [  9.49406E-02 0.00077 -1.39484E-03 0.00566  5.06551E-04 0.02362  8.46114E-02 0.00298 ];
INF_S3                    (idx, [1:   8]) = [  9.05804E-03 0.00573 -1.61223E-03 0.00363  1.81933E-04 0.03182  2.57982E-02 0.00790 ];
INF_S4                    (idx, [1:   8]) = [ -8.86825E-03 0.00467 -5.22091E-04 0.01527  2.46832E-06 1.00000 -5.23723E-03 0.02963 ];
INF_S5                    (idx, [1:   8]) = [  5.05213E-04 0.07100  2.25672E-05 0.22331 -6.50220E-05 0.08963  4.39116E-03 0.03727 ];
INF_S6                    (idx, [1:   8]) = [  5.20656E-03 0.00759 -1.24042E-04 0.04637 -9.29214E-05 0.07050 -1.23696E-02 0.01250 ];
INF_S7                    (idx, [1:   8]) = [  8.59524E-04 0.04433 -1.49383E-04 0.04080 -9.99698E-05 0.05103 -6.30299E-04 0.22968 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18675E-01 0.00027  1.58212E-02 0.00088  2.29507E-03 0.00628  1.22068E+00 0.00056 ];
INF_SP1                   (idx, [1:   8]) = [  2.31439E-01 0.00050  4.57101E-03 0.00225  9.22498E-04 0.01257  3.30826E-01 0.00124 ];
INF_SP2                   (idx, [1:   8]) = [  9.49409E-02 0.00077 -1.39484E-03 0.00566  5.06551E-04 0.02362  8.46114E-02 0.00298 ];
INF_SP3                   (idx, [1:   8]) = [  9.05820E-03 0.00573 -1.61223E-03 0.00363  1.81933E-04 0.03182  2.57982E-02 0.00790 ];
INF_SP4                   (idx, [1:   8]) = [ -8.86861E-03 0.00466 -5.22091E-04 0.01527  2.46832E-06 1.00000 -5.23723E-03 0.02963 ];
INF_SP5                   (idx, [1:   8]) = [  5.04718E-04 0.07112  2.25672E-05 0.22331 -6.50220E-05 0.08963  4.39116E-03 0.03727 ];
INF_SP6                   (idx, [1:   8]) = [  5.20659E-03 0.00758 -1.24042E-04 0.04637 -9.29214E-05 0.07050 -1.23696E-02 0.01250 ];
INF_SP7                   (idx, [1:   8]) = [  8.59450E-04 0.04450 -1.49383E-04 0.04080 -9.99698E-05 0.05103 -6.30299E-04 0.22968 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42160E-01 0.00134  7.87286E-01 0.00575 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.41567E-01 0.00169  7.87707E-01 0.00840 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.42129E-01 0.00190  7.85040E-01 0.01179 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.42821E-01 0.00215  7.93450E-01 0.01057 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37656E+00 0.00134  4.23730E-01 0.00572 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37998E+00 0.00169  4.23896E-01 0.00851 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37680E+00 0.00190  4.26041E-01 0.01192 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37291E+00 0.00215  4.21253E-01 0.01075 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.08126E-03 0.01617  2.16241E-04 0.09286  9.67851E-04 0.04562  6.48742E-04 0.05413  1.37145E-03 0.04156  2.30800E-03 0.02826  7.06620E-04 0.05239  6.05925E-04 0.05564  2.56431E-04 0.08624 ];
LAMBDA                    (idx, [1:  18]) = [  4.71813E-01 0.02617  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.6E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu Feb  8 00:47:51 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.3E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.69182E-02 0.00143  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73082E-01 4.0E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.34646E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98765E-01 7.5E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98260E-01 1.1E-05 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.34335E-01 0.00013  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.82485E+00 0.00050  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.68011E+01 0.00052  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.68011E+01 0.00052  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.65115E+00 0.00069  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.21629E-01 0.00155  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000518 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00104E+03 0.00143 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00104E+03 0.00143 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  3.53210E+02 ;
RUNNING_TIME              (idx, 1)        =  3.58457E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.18593E+00  2.33383E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  3.54261E+02  3.23368E+01  3.11457E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  1.34400E-01  6.78333E-03 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  5.81667E-02  1.51667E-03 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  3.58457E+02  3.42913E+03 ];
CPU_USAGE                 (idx, 1)        = 0.98536 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.98715E-01 0.00014 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.81618E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.10000E-06 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 93 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  4.08287E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.34604E+01 ;
TOT_SF_RATE               (idx, 1)        =  4.49443E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  5.94539E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  4.30659E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  3.48832E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.30296E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  1.26931E+04 ;
INGESTION_TOXICITY        (idx, 1)        =  1.42725E+04 ;
ACTINIDE_INH_TOX          (idx, 1)        =  3.36735E+03 ;
ACTINIDE_ING_TOX          (idx, 1)        =  2.43865E+03 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  9.32574E+03 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  1.18338E+04 ;
SR90_ACTIVITY             (idx, 1)        =  3.54497E+08 ;
TE132_ACTIVITY            (idx, 1)        =  2.69005E+11 ;
I131_ACTIVITY             (idx, 1)        =  1.35718E+11 ;
I132_ACTIVITY             (idx, 1)        =  2.70055E+11 ;
CS134_ACTIVITY            (idx, 1)        =  3.35959E+06 ;
CS137_ACTIVITY            (idx, 1)        =  3.78909E+08 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  4.48336E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.93754E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.59999E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  6.62156E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.91294E+09 0.00101  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 6 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  5.00000E-01  5.00187E-01 ];
BURN_DAYS                 (idx, [1:  2])  = [  1.48943E+01  4.46828E+00 ];
FIMA                      (idx, [1:  3])  = [  5.26208E-04  8.36279E+18  1.58842E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.96615E-01 0.00233 ];
U235_FISS                 (idx, [1:   4]) = [  6.07823E+12 0.00123  9.34420E-01 0.00037 ];
U238_FISS                 (idx, [1:   4]) = [  3.60526E+11 0.00603  5.54033E-02 0.00577 ];
PU239_FISS                (idx, [1:   4]) = [  6.49427E+10 0.01325  9.99469E-03 0.01335 ];
U235_CAPT                 (idx, [1:   4]) = [  1.45748E+12 0.00291  2.73050E-01 0.00245 ];
U238_CAPT                 (idx, [1:   4]) = [  2.99042E+12 0.00226  5.60080E-01 0.00123 ];
PU239_CAPT                (idx, [1:   4]) = [  3.74866E+10 0.01817  7.02209E-03 0.01806 ];
PU240_CAPT                (idx, [1:   4]) = [  1.29553E+09 0.09890  2.43002E-04 0.09875 ];
XE135_CAPT                (idx, [1:   4]) = [  2.63466E+11 0.00654  4.93801E-02 0.00649 ];
SM149_CAPT                (idx, [1:   4]) = [  4.64423E+10 0.01549  8.70111E-03 0.01544 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000518 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.49163E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000518 1.00149E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 450903 4.51360E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 549615 5.50132E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000518 1.00149E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -2.91038E-09 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.60008E+13 2.3E-05  1.60008E+13 2.3E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49520E+12 1.9E-06  6.49520E+12 1.9E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  5.32982E+12 0.00098  4.87016E+12 0.00104  4.59663E+11 0.00127 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.18250E+13 0.00044  1.13654E+13 0.00045  4.59663E+11 0.00127 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.18259E+13 0.00101  1.18259E+13 0.00101  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.97189E+14 0.00088  1.65112E+14 0.00097  3.32076E+14 0.00090 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.18250E+13 0.00044 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.16969E+14 0.00075 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27511E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27511E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.80854E+00 0.00080 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.54278E-01 0.00031 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.76601E-01 0.00088 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36238E+00 0.00082 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.35533E+00 0.00097 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.35533E+00 0.00097 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46348E+00 2.4E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02533E+02 1.9E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.35473E+00 0.00101  1.34586E+00 0.00098  9.47186E-03 0.01630 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.35531E+00 0.00044 ];
COL_KEFF                  (idx, [1:   2]) = [  1.35371E+00 0.00100 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.35531E+00 0.00044 ];
ABS_KINF                  (idx, [1:   2]) = [  1.35531E+00 0.00044 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.69577E+01 0.00039 ];
IMP_ALF                   (idx, [1:   2]) = [  1.69542E+01 0.00018 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.73393E-07 0.00673 ];
IMP_EALF                  (idx, [1:   2]) = [  8.68763E-07 0.00303 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.98280E-01 0.00575 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.97274E-01 0.00235 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.32295E-03 0.01244  1.64069E-04 0.06666  7.54311E-04 0.03202  4.89296E-04 0.03909  1.00152E-03 0.02769  1.69767E-03 0.02224  5.60507E-04 0.03728  4.88173E-04 0.03818  1.67407E-04 0.06623 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.57416E-01 0.01950  4.48801E-03 0.05969  2.46138E-02 0.01730  3.09578E-02 0.02736  1.23995E-01 0.01209  2.87788E-01 0.00571  5.11863E-01 0.02460  1.20974E+00 0.02654  1.28677E+00 0.05943 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  6.89378E-03 0.01648  2.14406E-04 0.09203  9.53579E-04 0.04290  6.35770E-04 0.05260  1.31488E-03 0.03821  2.19473E-03 0.03189  7.29894E-04 0.05171  6.20281E-04 0.05258  2.30246E-04 0.09735 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.55974E-01 0.02689  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.7E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.21199E-05 0.00212  1.21133E-05 0.00212  1.29193E-05 0.02270 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.64107E-05 0.00184  1.64017E-05 0.00185  1.74882E-05 0.02264 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  6.96995E-03 0.01687  2.32492E-04 0.09752  9.68934E-04 0.04252  6.67217E-04 0.05065  1.34904E-03 0.03668  2.18087E-03 0.03059  7.32552E-04 0.04995  6.14653E-04 0.05617  2.24194E-04 0.09006 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.44749E-01 0.02850  1.24667E-02 0.0E+00  2.82917E-02 1.7E-09  4.25244E-02 7.9E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.4E-09  3.55460E+00 5.0E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.20988E-05 0.00451  1.20911E-05 0.00453  9.82757E-06 0.05499 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.63825E-05 0.00441  1.63722E-05 0.00442  1.32991E-05 0.05483 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  7.08575E-03 0.04690  2.05527E-04 0.22417  1.06027E-03 0.11866  7.31740E-04 0.14573  1.21490E-03 0.11397  2.22821E-03 0.08638  8.18428E-04 0.13701  6.07910E-04 0.15636  2.18765E-04 0.27983 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  4.29549E-01 0.06868  1.24667E-02 4.6E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.0E-09  2.92467E-01 5.9E-09  6.66488E-01 5.3E-09  1.63478E+00 0.0E+00  3.55460E+00 4.7E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  7.10052E-03 0.04648  2.20442E-04 0.22359  1.05401E-03 0.11610  7.25507E-04 0.14276  1.22811E-03 0.11007  2.23711E-03 0.08279  8.08945E-04 0.13941  6.07738E-04 0.15586  2.18653E-04 0.25707 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.27018E-01 0.06769  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 4.0E-09  2.92467E-01 5.9E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 4.7E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -5.87702E+02 0.04687 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.21072E-05 0.00132 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.63940E-05 0.00090 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  7.05094E-03 0.00851 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.82756E+02 0.00856 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.09276E-07 0.00114 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92177E-06 0.00100  2.92204E-06 0.00101  2.87533E-06 0.01239 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.90150E-05 0.00129  1.90150E-05 0.00130  1.91528E-05 0.01662 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.77279E-01 0.00088  5.76081E-01 0.00088  8.77796E-01 0.02184 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.27530E+01 0.02683 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.67680E+01 0.00052  2.91978E+01 0.00067 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.30195E+04 0.00573  5.43458E+04 0.00358  1.14551E+05 0.00174  1.25943E+05 0.00182  1.18329E+05 0.00144  1.32938E+05 0.00141  9.03204E+04 0.00158  8.11342E+04 0.00188  6.16876E+04 0.00228  5.02814E+04 0.00196  4.31800E+04 0.00199  3.89190E+04 0.00171  3.59631E+04 0.00177  3.42044E+04 0.00208  3.31152E+04 0.00175  2.85619E+04 0.00290  2.82232E+04 0.00180  2.79675E+04 0.00170  2.72261E+04 0.00217  5.27892E+04 0.00163  5.02078E+04 0.00120  3.55672E+04 0.00248  2.27489E+04 0.00280  2.58399E+04 0.00193  2.40134E+04 0.00240  2.17219E+04 0.00286  3.47902E+04 0.00174  8.04063E+03 0.00496  1.01273E+04 0.00344  9.27882E+03 0.00336  5.30858E+03 0.00607  9.30788E+03 0.00368  6.31613E+03 0.00362  5.22829E+03 0.00448  9.86371E+02 0.00846  9.94425E+02 0.00877  9.96478E+02 0.00829  1.04030E+03 0.00862  1.03200E+03 0.00929  9.97923E+02 0.00888  1.05107E+03 0.00907  9.86338E+02 0.00713  1.87066E+03 0.00683  2.99470E+03 0.00501  3.78981E+03 0.00600  9.99265E+03 0.00289  1.03579E+04 0.00330  1.09000E+04 0.00269  6.87604E+03 0.00410  4.75028E+03 0.00491  3.56252E+03 0.00357  4.00291E+03 0.00397  7.12509E+03 0.00360  8.86509E+03 0.00321  1.58495E+04 0.00252  2.14347E+04 0.00255  2.85528E+04 0.00227  1.68523E+04 0.00306  1.15149E+04 0.00318  8.04375E+03 0.00352  7.05317E+03 0.00358  6.71130E+03 0.00369  5.53661E+03 0.00381  3.68798E+03 0.00394  3.35026E+03 0.00339  2.92998E+03 0.00483  2.45615E+03 0.00552  1.88328E+03 0.00507  1.23364E+03 0.00622  4.25024E+02 0.01069 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.35345E+00 0.00106 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.39902E+14 0.00105  5.73404E+13 0.00091 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.45926E-01 0.00027  1.34061E+00 0.00044 ];
INF_CAPT                  (idx, [1:   4]) = [  7.57919E-03 0.00117  3.48223E-02 0.00061 ];
INF_ABS                   (idx, [1:   4]) = [  1.14145E-02 0.00086  1.18729E-01 0.00072 ];
INF_FISS                  (idx, [1:   4]) = [  3.83534E-03 0.00094  8.39066E-02 0.00078 ];
INF_NSF                   (idx, [1:   4]) = [  9.68989E-03 0.00092  2.04847E-01 0.00078 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52648E+00 7.9E-05  2.44137E+00 2.1E-06 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03086E+02 5.5E-06  2.02339E+02 3.5E-07 ];
INF_INVV                  (idx, [1:   4]) = [  5.49125E-08 0.00087  2.26071E-06 0.00053 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34504E-01 0.00029  1.22164E+00 0.00049 ];
INF_SCATT1                (idx, [1:   4]) = [  2.36072E-01 0.00037  3.31185E-01 0.00088 ];
INF_SCATT2                (idx, [1:   4]) = [  9.35549E-02 0.00060  8.48148E-02 0.00252 ];
INF_SCATT3                (idx, [1:   4]) = [  7.35293E-03 0.00712  2.57765E-02 0.00837 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.43483E-03 0.00529 -5.20899E-03 0.03183 ];
INF_SCATT5                (idx, [1:   4]) = [  4.77401E-04 0.09240  4.68080E-03 0.02954 ];
INF_SCATT6                (idx, [1:   4]) = [  5.05131E-03 0.00775 -1.21236E-02 0.01119 ];
INF_SCATT7                (idx, [1:   4]) = [  8.39875E-04 0.04175 -5.21822E-04 0.23347 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34544E-01 0.00029  1.22164E+00 0.00049 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.36074E-01 0.00037  3.31185E-01 0.00088 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.35544E-02 0.00060  8.48148E-02 0.00252 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.35280E-03 0.00712  2.57765E-02 0.00837 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.43476E-03 0.00531 -5.20899E-03 0.03183 ];
INF_SCATTP5               (idx, [1:   4]) = [  4.77620E-04 0.09245  4.68080E-03 0.02954 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.05164E-03 0.00775 -1.21236E-02 0.01119 ];
INF_SCATTP7               (idx, [1:   4]) = [  8.39543E-04 0.04171 -5.21822E-04 0.23347 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30760E-01 0.00063  8.86966E-01 0.00064 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44452E+00 0.00063  3.75817E-01 0.00064 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13744E-02 0.00087  1.18729E-01 0.00072 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72241E-02 0.00039  1.21284E-01 0.00113 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18702E-01 0.00028  1.58028E-02 0.00089  2.31701E-03 0.00941  1.21933E+00 0.00049 ];
INF_S1                    (idx, [1:   8]) = [  2.31506E-01 0.00038  4.56670E-03 0.00147  9.22017E-04 0.01640  3.30263E-01 0.00088 ];
INF_S2                    (idx, [1:   8]) = [  9.49369E-02 0.00058 -1.38203E-03 0.00450  4.94760E-04 0.02608  8.43200E-02 0.00250 ];
INF_S3                    (idx, [1:   8]) = [  8.96600E-03 0.00568 -1.61307E-03 0.00381  1.73010E-04 0.04495  2.56035E-02 0.00842 ];
INF_S4                    (idx, [1:   8]) = [ -8.90352E-03 0.00541 -5.31305E-04 0.00888 -2.68815E-06 1.00000 -5.20630E-03 0.03184 ];
INF_S5                    (idx, [1:   8]) = [  4.60104E-04 0.09377  1.72969E-05 0.38568 -6.91423E-05 0.12108  4.74994E-03 0.02989 ];
INF_S6                    (idx, [1:   8]) = [  5.17035E-03 0.00750 -1.19043E-04 0.04213 -8.36451E-05 0.08644 -1.20399E-02 0.01135 ];
INF_S7                    (idx, [1:   8]) = [  9.82727E-04 0.03664 -1.42852E-04 0.02867 -7.85381E-05 0.07299 -4.43284E-04 0.27321 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18742E-01 0.00028  1.58028E-02 0.00089  2.31701E-03 0.00941  1.21933E+00 0.00049 ];
INF_SP1                   (idx, [1:   8]) = [  2.31507E-01 0.00038  4.56670E-03 0.00147  9.22017E-04 0.01640  3.30263E-01 0.00088 ];
INF_SP2                   (idx, [1:   8]) = [  9.49364E-02 0.00059 -1.38203E-03 0.00450  4.94760E-04 0.02608  8.43200E-02 0.00250 ];
INF_SP3                   (idx, [1:   8]) = [  8.96587E-03 0.00569 -1.61307E-03 0.00381  1.73010E-04 0.04495  2.56035E-02 0.00842 ];
INF_SP4                   (idx, [1:   8]) = [ -8.90346E-03 0.00543 -5.31305E-04 0.00888 -2.68815E-06 1.00000 -5.20630E-03 0.03184 ];
INF_SP5                   (idx, [1:   8]) = [  4.60323E-04 0.09379  1.72969E-05 0.38568 -6.91423E-05 0.12108  4.74994E-03 0.02989 ];
INF_SP6                   (idx, [1:   8]) = [  5.17068E-03 0.00750 -1.19043E-04 0.04213 -8.36451E-05 0.08644 -1.20399E-02 0.01135 ];
INF_SP7                   (idx, [1:   8]) = [  9.82395E-04 0.03662 -1.42852E-04 0.02867 -7.85381E-05 0.07299 -4.43284E-04 0.27321 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42840E-01 0.00123  7.91864E-01 0.00613 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.42949E-01 0.00207  7.78338E-01 0.00869 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.42544E-01 0.00210  8.08984E-01 0.00903 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.43060E-01 0.00128  7.91373E-01 0.00938 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37270E+00 0.00124  4.21332E-01 0.00620 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37217E+00 0.00207  4.29044E-01 0.00874 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37447E+00 0.00210  4.12857E-01 0.00916 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37146E+00 0.00128  4.22095E-01 0.00934 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  6.89378E-03 0.01648  2.14406E-04 0.09203  9.53579E-04 0.04290  6.35770E-04 0.05260  1.31488E-03 0.03821  2.19473E-03 0.03189  7.29894E-04 0.05171  6.20281E-04 0.05258  2.30246E-04 0.09735 ];
LAMBDA                    (idx, [1:  18]) = [  4.55974E-01 0.02689  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.7E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];


% Increase counter:

if (exist('idx', 'var'));
  idx = idx + 1;
else;
  idx = 1;
end;

% Version, title and date:

VERSION                   (idx, [1: 14])  = 'Serpent 2.1.32' ;
COMPILE_DATE              (idx, [1: 20])  = 'Feb  6 2024 17:46:48' ;
DEBUG                     (idx, 1)        = 0 ;
TITLE                     (idx, [1:  8])  = 'Untitled' ;
CONFIDENTIAL_DATA         (idx, 1)        = 0 ;
INPUT_FILE_NAME           (idx, [1: 18])  = 'TMI_4rings_jeff311' ;
WORKING_DIRECTORY         (idx, [1: 35])  = '/home/p117902/Serpent2/Linux_x86_64' ;
HOSTNAME                  (idx, [1: 28])  = 'doppler.recherche.polymtl.ca' ;
CPU_TYPE                  (idx, [1: 47])  = 'Intel(R) Xeon(R) CPU           X5650  @ 2.67GHz' ;
CPU_MHZ                   (idx, 1)        = 31.0 ;
START_DATE                (idx, [1: 24])  = 'Wed Feb  7 18:49:24 2024' ;
COMPLETE_DATE             (idx, [1: 24])  = 'Thu Feb  8 01:50:06 2024' ;

% Run parameters:

POP                       (idx, 1)        = 2000 ;
CYCLES                    (idx, 1)        = 500 ;
SKIP                      (idx, 1)        = 20 ;
BATCH_INTERVAL            (idx, 1)        = 1 ;
SRC_NORM_MODE             (idx, 1)        = 2 ;
SEED                      (idx, 1)        = 1707349764305 ;
UFS_MODE                  (idx, 1)        = 0 ;
UFS_ORDER                 (idx, 1)        = 1.00000;
NEUTRON_TRANSPORT_MODE    (idx, 1)        = 1 ;
PHOTON_TRANSPORT_MODE     (idx, 1)        = 0 ;
GROUP_CONSTANT_GENERATION (idx, 1)        = 1 ;
B1_CALCULATION            (idx, [1:  3])  = [ 0 0 0 ];
B1_BURNUP_CORRECTION      (idx, 1)        = 0 ;

CRIT_SPEC_MODE            (idx, 1)        = 0 ;
IMPLICIT_REACTION_RATES   (idx, 1)        = 1 ;

% Optimization:

OPTIMIZATION_MODE         (idx, 1)        = 4 ;
RECONSTRUCT_MICROXS       (idx, 1)        = 1 ;
RECONSTRUCT_MACROXS       (idx, 1)        = 1 ;
DOUBLE_INDEXING           (idx, 1)        = 0 ;
MG_MAJORANT_MODE          (idx, 1)        = 0 ;
SPECTRUM_COLLAPSE         (idx, 1)        = 0 ;

% Parallelization:

MPI_TASKS                 (idx, 1)        = 1 ;
OMP_THREADS               (idx, 1)        = 1 ;
MPI_REPRODUCIBILITY       (idx, 1)        = 0 ;
OMP_REPRODUCIBILITY       (idx, 1)        = 1 ;
SHARE_BUF_ARRAY           (idx, 1)        = 1 ;
SHARE_RES2_ARRAY          (idx, 1)        = 1 ;
OMP_SHARED_QUEUE_LIM      (idx, 1)        = 0 ;

% File paths:

XS_DATA_FILE_PATH         (idx, [1: 31])  = '../xs/jeff311/sss_jeff311u.data' ;
DECAY_DATA_FILE_PATH      (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.dec' ;
SFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
NFY_DATA_FILE_PATH        (idx, [1: 29])  = '../xs/jeff311/sss_jeff311.nfy' ;
BRA_DATA_FILE_PATH        (idx, [1:  3])  = 'N/A' ;

% Collision and reaction sampling (neutrons/photons):

MIN_MACROXS               (idx, [1:   4]) = [  5.00000E-02 1.1E-09  0.00000E+00 0.0E+00 ];
DT_THRESH                 (idx, [1:  2])  = [  9.00000E-01  9.00000E-01 ];
ST_FRAC                   (idx, [1:   4]) = [  2.68260E-02 0.00142  0.00000E+00 0.0E+00 ];
DT_FRAC                   (idx, [1:   4]) = [  9.73174E-01 3.9E-05  0.00000E+00 0.0E+00 ];
DT_EFF                    (idx, [1:   4]) = [  7.34886E-01 0.00013  0.00000E+00 0.0E+00 ];
REA_SAMPLING_EFF          (idx, [1:   4]) = [  9.98789E-01 6.8E-06  0.00000E+00 0.0E+00 ];
REA_SAMPLING_FAIL         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_SAMPLING_EFF          (idx, [1:   2]) = [  9.98293E-01 9.5E-06 ];
TOT_COL_EFF               (idx, [1:   4]) = [  7.34586E-01 0.00013  0.00000E+00 0.0E+00 ];
AVG_TRACKING_LOOPS        (idx, [1:   8]) = [  2.82287E+00 0.00049  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TMS_FAIL_STAT             (idx, [1:   6]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
AVG_TRACKS                (idx, [1:   4]) = [  2.67878E+01 0.00054  0.00000E+00 0.0E+00 ];
AVG_REAL_COL              (idx, [1:   4]) = [  2.67878E+01 0.00054  0.00000E+00 0.0E+00 ];
AVG_VIRT_COL              (idx, [1:   4]) = [  9.63488E+00 0.00073  0.00000E+00 0.0E+00 ];
AVG_SURF_CROSS            (idx, [1:   4]) = [  9.17922E-01 0.00152  0.00000E+00 0.0E+00 ];
LOST_PARTICLES            (idx, 1)        = 0 ;

% Run statistics:

CYCLE_IDX                 (idx, 1)        = 500 ;
SIMULATED_HISTORIES       (idx, 1)        = 1000406 ;
MEAN_POP_SIZE             (idx, [1:  2])  = [  2.00081E+03 0.00128 ];
MEAN_POP_WGT              (idx, [1:  2])  = [  2.00081E+03 0.00128 ];
SIMULATION_COMPLETED      (idx, 1)        = 1 ;

% Running times:

TOT_CPU_TIME              (idx, 1)        =  4.15222E+02 ;
RUNNING_TIME              (idx, 1)        =  4.20698E+02 ;
INIT_TIME                 (idx, [1:  2])  = [  1.81308E+00  1.81308E+00 ];
PROCESS_TIME              (idx, [1:  2])  = [  2.62062E+00  2.24383E-01 ];
TRANSPORT_CYCLE_TIME      (idx, [1:  3])  = [  4.16042E+02  3.11035E+01  3.06775E+01 ];
BURNUP_CYCLE_TIME         (idx, [1:  2])  = [  1.46483E-01  5.88333E-03 ];
BATEMAN_SOLUTION_TIME     (idx, [1:  2])  = [  6.87667E-02  1.05000E-03 ];
MPI_OVERHEAD_TIME         (idx, [1:  2])  = [  0.00000E+00  0.00000E+00 ];
ESTIMATED_RUNNING_TIME    (idx, [1:  2])  = [  4.20689E+02  3.46748E+03 ];
CPU_USAGE                 (idx, 1)        = 0.98698 ;
TRANSPORT_CPU_USAGE       (idx, [1:   2]) = [  9.99034E-01 0.00011 ];
OMP_PARALLEL_FRAC         (idx, 1)        =  9.83899E-01 ;

% Memory usage:

AVAIL_MEM                 (idx, 1)        = 161016.26 ;
ALLOC_MEMSIZE             (idx, 1)        = 17616.44;
MEMSIZE                   (idx, 1)        = 17562.06;
XS_MEMSIZE                (idx, 1)        = 17293.81;
MAT_MEMSIZE               (idx, 1)        = 253.88;
RES_MEMSIZE               (idx, 1)        = 0.92;
IFC_MEMSIZE               (idx, 1)        = 0.00;
MISC_MEMSIZE              (idx, 1)        = 13.45;
UNKNOWN_MEMSIZE           (idx, 1)        = 0.00;
UNUSED_MEMSIZE            (idx, 1)        = 54.38;

% Geometry parameters:

TOT_CELLS                 (idx, 1)        = 8 ;
UNION_CELLS               (idx, 1)        = 0 ;

% Neutron energy grid:

NEUTRON_ERG_TOL           (idx, 1)        =  0.00000E+00 ;
NEUTRON_ERG_NE            (idx, 1)        = 940251 ;
NEUTRON_EMIN              (idx, 1)        =  1.00000E-11 ;
NEUTRON_EMAX              (idx, 1)        =  2.00000E+01 ;

% Unresolved resonance probability table sampling:

URES_DILU_CUT             (idx, 1)        =  1.00000E-09 ;
URES_EMIN                 (idx, 1)        =  1.10000E-06 ;
URES_EMAX                 (idx, 1)        =  1.00000E+00 ;
URES_AVAIL                (idx, 1)        = 162 ;
URES_USED                 (idx, 1)        = 96 ;

% Nuclides and reaction channels:

TOT_NUCLIDES              (idx, 1)        = 1458 ;
TOT_TRANSPORT_NUCLIDES    (idx, 1)        = 328 ;
TOT_DOSIMETRY_NUCLIDES    (idx, 1)        = 0 ;
TOT_DECAY_NUCLIDES        (idx, 1)        = 1130 ;
TOT_PHOTON_NUCLIDES       (idx, 1)        = 0 ;
TOT_REA_CHANNELS          (idx, 1)        = 7283 ;
TOT_TRANSMU_REA           (idx, 1)        = 2340 ;

% Neutron physics options:

USE_DELNU                 (idx, 1)        = 1 ;
USE_URES                  (idx, 1)        = 1 ;
USE_DBRC                  (idx, 1)        = 0 ;
IMPL_CAPT                 (idx, 1)        = 0 ;
IMPL_NXN                  (idx, 1)        = 1 ;
IMPL_FISS                 (idx, 1)        = 0 ;
DOPPLER_PREPROCESSOR      (idx, 1)        = 0 ;
TMS_MODE                  (idx, 1)        = 2 ;
SAMPLE_FISS               (idx, 1)        = 1 ;
SAMPLE_CAPT               (idx, 1)        = 1 ;
SAMPLE_SCATT              (idx, 1)        = 1 ;

% Energy deposition:

EDEP_MODE                 (idx, 1)        = 0 ;
EDEP_DELAYED              (idx, 1)        = 1 ;
EDEP_KEFF_CORR            (idx, 1)        = 1 ;
EDEP_LOCAL_EGD            (idx, 1)        = 0 ;
EDEP_COMP                 (idx, [1:  9])  = [ 0 0 0 0 0 0 0 0 0 ];
EDEP_CAPT_E               (idx, 1)        =  0.00000E+00 ;

% Radioactivity data:

TOT_ACTIVITY              (idx, 1)        =  4.13802E+13 ;
TOT_DECAY_HEAT            (idx, 1)        =  1.35162E+01 ;
TOT_SF_RATE               (idx, 1)        =  5.11896E-02 ;
ACTINIDE_ACTIVITY         (idx, 1)        =  6.00529E+12 ;
ACTINIDE_DECAY_HEAT       (idx, 1)        =  4.34925E-01 ;
FISSION_PRODUCT_ACTIVITY  (idx, 1)        =  3.53748E+13 ;
FISSION_PRODUCT_DECAY_HEAT(idx, 1)        =  1.30812E+01 ;
INHALATION_TOXICITY       (idx, 1)        =  1.50323E+04 ;
INGESTION_TOXICITY        (idx, 1)        =  1.55665E+04 ;
ACTINIDE_INH_TOX          (idx, 1)        =  3.63818E+03 ;
ACTINIDE_ING_TOX          (idx, 1)        =  2.47851E+03 ;
FISSION_PRODUCT_INH_TOX   (idx, 1)        =  1.13941E+04 ;
FISSION_PRODUCT_ING_TOX   (idx, 1)        =  1.30880E+04 ;
SR90_ACTIVITY             (idx, 1)        =  5.30520E+08 ;
TE132_ACTIVITY            (idx, 1)        =  2.78243E+11 ;
I131_ACTIVITY             (idx, 1)        =  1.61512E+11 ;
I132_ACTIVITY             (idx, 1)        =  2.79628E+11 ;
CS134_ACTIVITY            (idx, 1)        =  1.03049E+07 ;
CS137_ACTIVITY            (idx, 1)        =  5.68339E+08 ;
PHOTON_DECAY_SOURCE       (idx, 1)        =  4.53463E+13 ;
NEUTRON_DECAY_SOURCE      (idx, 1)        =  8.90816E+10 ;
ALPHA_DECAY_SOURCE        (idx, 1)        =  1.77448E+07 ;
ELECTRON_DECAY_SOURCE     (idx, 1)        =  6.73235E+13 ;

% Normalization coefficient:

NORM_COEF                 (idx, [1:   4]) = [  5.93206E+09 0.00099  0.00000E+00 0.0E+00 ];

% Parameters for burnup calculation:

BURN_MATERIALS            (idx, 1)        = 4 ;
BURN_MODE                 (idx, 1)        = 2 ;
BURN_STEP                 (idx, 1)        = 7 ;
BURN_RANDOMIZE_DATA       (idx, [1:  3])  = [ 0 0 0 ];
BURNUP                    (idx, [1:  2])  = [  7.50000E-01  7.50290E-01 ];
BURN_DAYS                 (idx, [1:  2])  = [  2.23414E+01  7.44713E+00 ];
FIMA                      (idx, [1:  3])  = [  7.89262E-04  1.25434E+19  1.58800E+22 ];

% Analog reaction rate estimators:

CONVERSION_RATIO          (idx, [1:   2]) = [  3.95961E-01 0.00238 ];
U235_FISS                 (idx, [1:   4]) = [  6.04032E+12 0.00113  9.28309E-01 0.00037 ];
U238_FISS                 (idx, [1:   4]) = [  3.60712E+11 0.00580  5.53981E-02 0.00545 ];
PU239_FISS                (idx, [1:   4]) = [  1.04616E+11 0.00993  1.60758E-02 0.00985 ];
PU240_FISS                (idx, [1:   4]) = [  1.13254E+07 1.00000  1.75439E-06 1.00000 ];
PU241_FISS                (idx, [1:   4]) = [  3.66073E+07 0.57629  5.54955E-06 0.57624 ];
U235_CAPT                 (idx, [1:   4]) = [  1.44828E+12 0.00277  2.69552E-01 0.00243 ];
U238_CAPT                 (idx, [1:   4]) = [  2.98926E+12 0.00230  5.56126E-01 0.00136 ];
PU239_CAPT                (idx, [1:   4]) = [  6.01486E+10 0.01459  1.11966E-02 0.01453 ];
PU240_CAPT                (idx, [1:   4]) = [  3.43783E+09 0.05941  6.41634E-04 0.05971 ];
PU241_CAPT                (idx, [1:   4]) = [  1.24029E+07 1.00000  2.33372E-06 1.00000 ];
XE135_CAPT                (idx, [1:   4]) = [  2.60961E+11 0.00679  4.85693E-02 0.00666 ];
SM149_CAPT                (idx, [1:   4]) = [  5.79521E+10 0.01455  1.07957E-02 0.01463 ];

% Neutron balance (particles/weight):

BALA_SRC_NEUTRON_SRC     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_FISS    (idx, [1:  2])  = [ 1000406 1.00000E+06 ];
BALA_SRC_NEUTRON_NXN     (idx, [1:  2])  = [ 0 1.44345E+03 ];
BALA_SRC_NEUTRON_VR      (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_SRC_NEUTRON_TOT     (idx, [1:  2])  = [ 1000406 1.00144E+06 ];

BALA_LOSS_NEUTRON_CAPT    (idx, [1:  2])  = [ 452414 4.52915E+05 ];
BALA_LOSS_NEUTRON_FISS    (idx, [1:  2])  = [ 547992 5.48528E+05 ];
BALA_LOSS_NEUTRON_LEAK    (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_CUT     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_ERR     (idx, [1:  2])  = [ 0 0.00000E+00 ];
BALA_LOSS_NEUTRON_TOT     (idx, [1:  2])  = [ 1000406 1.00144E+06 ];

BALA_NEUTRON_DIFF         (idx, [1:  2])  = [ 0 -8.14907E-10 ];

% Normalized total reaction rates (neutrons):

TOT_POWER                 (idx, [1:   6]) = [  2.10765E+02 5.1E-09  2.10765E+02 5.1E-09  0.00000E+00 0.0E+00 ];
TOT_POWDENS               (idx, [1:   6]) = [  3.35700E-02 4.5E-09  3.35700E-02 4.5E-09  0.00000E+00 0.0E+00 ];
TOT_GENRATE               (idx, [1:   6]) = [  1.60158E+13 2.4E-05  1.60158E+13 2.4E-05  0.00000E+00 0.0E+00 ];
TOT_FISSRATE              (idx, [1:   6]) = [  6.49407E+12 2.1E-06  6.49407E+12 2.1E-06  0.00000E+00 0.0E+00 ];
TOT_CAPTRATE              (idx, [1:   6]) = [  5.37656E+12 0.00092  4.91701E+12 0.00099  4.59548E+11 0.00119 ];
TOT_ABSRATE               (idx, [1:   6]) = [  1.18706E+13 0.00042  1.14111E+13 0.00043  4.59548E+11 0.00119 ];
TOT_SRCRATE               (idx, [1:   6]) = [  1.18641E+13 0.00099  1.18641E+13 0.00099  0.00000E+00 0.0E+00 ];
TOT_FLUX                  (idx, [1:   6]) = [  4.98484E+14 0.00087  1.65644E+14 0.00095  3.32839E+14 0.00089 ];
TOT_PHOTON_PRODRATE       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
TOT_LEAKRATE              (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
ALBEDO_LEAKRATE           (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_LOSSRATE              (idx, [1:   2]) = [  1.18706E+13 0.00042 ];
TOT_CUTRATE               (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
TOT_RR                    (idx, [1:   2]) = [  3.17814E+14 0.00071 ];
INI_FMASS                 (idx, 1)        =  6.27837E-03 ;
TOT_FMASS                 (idx, 1)        =  6.27348E-03 ;
INI_BURN_FMASS            (idx, 1)        =  6.27837E-03 ;
TOT_BURN_FMASS            (idx, 1)        =  6.27348E-03 ;

% Six-factor formula:

SIX_FF_ETA                (idx, [1:   2]) = [  1.80336E+00 0.00076 ];
SIX_FF_F                  (idx, [1:   2]) = [  9.54599E-01 0.00029 ];
SIX_FF_P                  (idx, [1:   2]) = [  5.76798E-01 0.00090 ];
SIX_FF_EPSILON            (idx, [1:   2]) = [  1.36285E+00 0.00084 ];
SIX_FF_LF                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_LT                 (idx, [1:   2]) = [  1.00000E+00 0.0E+00 ];
SIX_FF_KINF               (idx, [1:   2]) = [  1.35276E+00 0.00088 ];
SIX_FF_KEFF               (idx, [1:   2]) = [  1.35276E+00 0.00088 ];

% Fission neutron and energy production:

NUBAR                     (idx, [1:   2]) = [  2.46622E+00 2.6E-05 ];
FISSE                     (idx, [1:   2]) = [  2.02568E+02 2.1E-06 ];

% Criticality eigenvalues:

ANA_KEFF                  (idx, [1:   6]) = [  1.35274E+00 0.00092  1.34327E+00 0.00089  9.49092E-03 0.01599 ];
IMP_KEFF                  (idx, [1:   2]) = [  1.35141E+00 0.00042 ];
COL_KEFF                  (idx, [1:   2]) = [  1.35059E+00 0.00098 ];
ABS_KEFF                  (idx, [1:   2]) = [  1.35141E+00 0.00042 ];
ABS_KINF                  (idx, [1:   2]) = [  1.35141E+00 0.00042 ];
GEOM_ALBEDO               (idx, [1:   6]) = [  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00  1.00000E+00 0.0E+00 ];

% ALF (Average lethargy of neutrons causing fission):
% Based on E0 = 2.000000E+01 MeV

ANA_ALF                   (idx, [1:   2]) = [  1.69534E+01 0.00040 ];
IMP_ALF                   (idx, [1:   2]) = [  1.69482E+01 0.00019 ];

% EALF (Energy corresponding to average lethargy of neutrons causing fission):

ANA_EALF                  (idx, [1:   2]) = [  8.77664E-07 0.00688 ];
IMP_EALF                  (idx, [1:   2]) = [  8.74325E-07 0.00326 ];

% AFGE (Average energy of neutrons causing fission):

ANA_AFGE                  (idx, [1:   2]) = [  1.97716E-01 0.00568 ];
IMP_AFGE                  (idx, [1:   2]) = [  1.98124E-01 0.00251 ];

% Forward-weighted delayed neutron parameters:

PRECURSOR_GROUPS          (idx, 1)        = 8 ;
FWD_ANA_BETA_ZERO         (idx, [1:  18]) = [  5.44415E-03 0.01237  1.45755E-04 0.07533  8.13544E-04 0.03013  4.77421E-04 0.03968  1.03235E-03 0.02609  1.73634E-03 0.02122  5.60371E-04 0.03757  4.76287E-04 0.03973  2.02082E-04 0.05565 ];
FWD_ANA_LAMBDA            (idx, [1:  18]) = [  4.73018E-01 0.01819  3.86468E-03 0.06679  2.51230E-02 0.01590  3.10428E-02 0.02723  1.23995E-01 0.01209  2.85448E-01 0.00702  5.09197E-01 0.02488  1.20320E+00 0.02681  1.59246E+00 0.04969 ];

% Beta-eff using Meulekamp's method:

ADJ_MEULEKAMP_BETA_EFF    (idx, [1:  18]) = [  7.07783E-03 0.01686  1.93251E-04 0.09911  1.05645E-03 0.04354  6.29916E-04 0.05654  1.33656E-03 0.03645  2.29965E-03 0.02864  6.94495E-04 0.05078  6.16753E-04 0.05472  2.50759E-04 0.07982 ];
ADJ_MEULEKAMP_LAMBDA      (idx, [1:  18]) = [  4.63490E-01 0.02614  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.7E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];

% Adjoint weighted time constants using Nauchi's method:

IFP_CHAIN_LENGTH          (idx, 1)        = 15 ;
ADJ_NAUCHI_GEN_TIME       (idx, [1:   6]) = [  1.20761E-05 0.00216  1.20690E-05 0.00217  1.29008E-05 0.02140 ];
ADJ_NAUCHI_LIFETIME       (idx, [1:   6]) = [  1.63283E-05 0.00191  1.63187E-05 0.00192  1.74413E-05 0.02135 ];
ADJ_NAUCHI_BETA_EFF       (idx, [1:  18]) = [  7.00663E-03 0.01654  1.89349E-04 0.09863  1.04259E-03 0.04297  6.04342E-04 0.05504  1.34152E-03 0.03519  2.28657E-03 0.02868  6.98906E-04 0.05265  6.06665E-04 0.05529  2.36694E-04 0.08720 ];
ADJ_NAUCHI_LAMBDA         (idx, [1:  18]) = [  4.55699E-01 0.02740  1.24667E-02 0.0E+00  2.82917E-02 1.7E-09  4.25244E-02 8.1E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 5.5E-09  3.55460E+00 4.9E-09 ];

% Adjoint weighted time constants using IFP:

ADJ_IFP_GEN_TIME          (idx, [1:   6]) = [  1.20945E-05 0.00449  1.20881E-05 0.00450  9.60199E-06 0.05036 ];
ADJ_IFP_LIFETIME          (idx, [1:   6]) = [  1.63529E-05 0.00436  1.63443E-05 0.00437  1.29875E-05 0.05035 ];
ADJ_IFP_IMP_BETA_EFF      (idx, [1:  18]) = [  6.88333E-03 0.05263  1.86552E-04 0.28676  1.04341E-03 0.13074  5.96767E-04 0.15251  1.16366E-03 0.12773  2.13822E-03 0.08214  7.92046E-04 0.15859  6.29633E-04 0.15988  3.33045E-04 0.23245 ];
ADJ_IFP_IMP_LAMBDA        (idx, [1:  18]) = [  5.02382E-01 0.06715  1.24667E-02 2.7E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 3.7E-09  2.92467E-01 6.0E-09  6.66488E-01 5.3E-09  1.63478E+00 0.0E+00  3.55460E+00 6.0E-09 ];
ADJ_IFP_ANA_BETA_EFF      (idx, [1:  18]) = [  6.80607E-03 0.05006  1.86563E-04 0.27526  1.00846E-03 0.12054  6.27536E-04 0.14962  1.15751E-03 0.11941  2.07325E-03 0.08197  8.05923E-04 0.15057  5.96098E-04 0.15618  3.50724E-04 0.23584 ];
ADJ_IFP_ANA_LAMBDA        (idx, [1:  18]) = [  4.98075E-01 0.06787  1.24667E-02 3.8E-09  2.82917E-02 0.0E+00  4.25244E-02 0.0E+00  1.33042E-01 3.7E-09  2.92467E-01 6.0E-09  6.66488E-01 5.1E-09  1.63478E+00 0.0E+00  3.55460E+00 6.0E-09 ];
ADJ_IFP_ROSSI_ALPHA       (idx, [1:   2]) = [ -5.75255E+02 0.05234 ];

% Adjoint weighted time constants using perturbation technique:

ADJ_PERT_GEN_TIME         (idx, [1:   2]) = [  1.20813E-05 0.00123 ];
ADJ_PERT_LIFETIME         (idx, [1:   2]) = [  1.63363E-05 0.00085 ];
ADJ_PERT_BETA_EFF         (idx, [1:   2]) = [  6.93023E-03 0.00893 ];
ADJ_PERT_ROSSI_ALPHA      (idx, [1:   2]) = [ -5.73923E+02 0.00894 ];

% Inverse neutron speed :

ANA_INV_SPD               (idx, [1:   2]) = [  3.08067E-07 0.00116 ];

% Analog slowing-down and thermal neutron lifetime (total/prompt/delayed):

ANA_SLOW_TIME             (idx, [1:   6]) = [  2.92324E-06 0.00098  2.92348E-06 0.00099  2.88208E-06 0.01271 ];
ANA_THERM_TIME            (idx, [1:   6]) = [  1.89036E-05 0.00121  1.89008E-05 0.00121  1.93002E-05 0.01528 ];
ANA_THERM_FRAC            (idx, [1:   6]) = [  5.77462E-01 0.00090  5.76217E-01 0.00092  8.85712E-01 0.02438 ];
ANA_DELAYED_EMTIME        (idx, [1:   2]) = [  1.17921E+01 0.02678 ];
ANA_MEAN_NCOL             (idx, [1:   4]) = [  2.67553E+01 0.00054  2.91770E+01 0.00068 ];

% Group constant generation:

GC_UNIVERSE_NAME          (idx, [1:  1])  = '0' ;

% Micro- and macro-group structures:

MICRO_NG                  (idx, 1)        = 70 ;
MICRO_E                   (idx, [1:  71]) = [  2.00000E+01  6.06550E+00  3.67900E+00  2.23100E+00  1.35300E+00  8.21000E-01  5.00000E-01  3.02500E-01  1.83000E-01  1.11000E-01  6.74300E-02  4.08500E-02  2.47800E-02  1.50300E-02  9.11800E-03  5.50000E-03  3.51910E-03  2.23945E-03  1.42510E-03  9.06898E-04  3.67262E-04  1.48728E-04  7.55014E-05  4.80520E-05  2.77000E-05  1.59680E-05  9.87700E-06  4.00000E-06  3.30000E-06  2.60000E-06  2.10000E-06  1.85500E-06  1.50000E-06  1.30000E-06  1.15000E-06  1.12300E-06  1.09700E-06  1.07100E-06  1.04500E-06  1.02000E-06  9.96000E-07  9.72000E-07  9.50000E-07  9.10000E-07  8.50000E-07  7.80000E-07  6.25000E-07  5.00000E-07  4.00000E-07  3.50000E-07  3.20000E-07  3.00000E-07  2.80000E-07  2.50000E-07  2.20000E-07  1.80000E-07  1.40000E-07  1.00000E-07  8.00000E-08  6.70000E-08  5.80000E-08  5.00000E-08  4.20000E-08  3.50000E-08  3.00000E-08  2.50000E-08  2.00000E-08  1.50000E-08  1.00000E-08  5.00000E-09  1.00000E-11 ];

MACRO_NG                  (idx, 1)        = 2 ;
MACRO_E                   (idx, [1:   3]) = [  1.00000E+37  6.25000E-07  0.00000E+00 ];

% Micro-group spectrum:

INF_MICRO_FLX             (idx, [1: 140]) = [  1.32481E+04 0.00940  5.39774E+04 0.00318  1.13980E+05 0.00181  1.25970E+05 0.00179  1.18485E+05 0.00134  1.32882E+05 0.00187  9.04283E+04 0.00140  8.14011E+04 0.00157  6.18314E+04 0.00180  5.03164E+04 0.00185  4.32568E+04 0.00198  3.90041E+04 0.00252  3.59889E+04 0.00176  3.41713E+04 0.00181  3.30797E+04 0.00205  2.86490E+04 0.00265  2.81810E+04 0.00228  2.78768E+04 0.00174  2.72120E+04 0.00267  5.27926E+04 0.00146  5.00694E+04 0.00155  3.56564E+04 0.00181  2.26808E+04 0.00214  2.58418E+04 0.00270  2.39576E+04 0.00188  2.17927E+04 0.00299  3.48173E+04 0.00183  8.02561E+03 0.00423  1.01288E+04 0.00402  9.21832E+03 0.00446  5.29984E+03 0.00513  9.26431E+03 0.00373  6.30084E+03 0.00528  5.26878E+03 0.00419  9.74725E+02 0.01017  9.78006E+02 0.00987  1.00338E+03 0.01035  1.03185E+03 0.00994  1.02676E+03 0.00733  1.01260E+03 0.01011  1.06490E+03 0.01103  9.93208E+02 0.00771  1.88667E+03 0.00691  2.98611E+03 0.00704  3.81513E+03 0.00469  9.98275E+03 0.00328  1.03987E+04 0.00362  1.08138E+04 0.00370  6.84005E+03 0.00416  4.70634E+03 0.00388  3.55362E+03 0.00618  3.94576E+03 0.00452  7.03340E+03 0.00386  8.77595E+03 0.00344  1.58130E+04 0.00243  2.14073E+04 0.00228  2.84071E+04 0.00233  1.67956E+04 0.00246  1.14400E+04 0.00308  8.01208E+03 0.00240  7.01042E+03 0.00330  6.70011E+03 0.00297  5.54005E+03 0.00310  3.65339E+03 0.00402  3.32184E+03 0.00431  2.90975E+03 0.00359  2.42185E+03 0.00515  1.88549E+03 0.00590  1.21301E+03 0.00384  4.17689E+02 0.00605 ];

% Integral parameters:

INF_KINF                  (idx, [1:   2]) = [  1.35093E+00 0.00087 ];

% Flux spectra in infinite geometry:

INF_FLX                   (idx, [1:   4]) = [  4.41289E+14 0.00089  5.72487E+13 0.00099 ];
INF_FISS_FLX              (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

INF_TOT                   (idx, [1:   4]) = [  5.46157E-01 0.00020  1.34259E+00 0.00046 ];
INF_CAPT                  (idx, [1:   4]) = [  7.59430E-03 0.00105  3.53964E-02 0.00077 ];
INF_ABS                   (idx, [1:   4]) = [  1.14197E-02 0.00078  1.19401E-01 0.00092 ];
INF_FISS                  (idx, [1:   4]) = [  3.82538E-03 0.00101  8.40049E-02 0.00099 ];
INF_NSF                   (idx, [1:   4]) = [  9.66976E-03 0.00101  2.05358E-01 0.00099 ];
INF_NUBAR                 (idx, [1:   4]) = [  2.52779E+00 7.6E-05  2.44459E+00 4.4E-06 ];
INF_KAPPA                 (idx, [1:   4]) = [  2.03099E+02 6.1E-06  2.02381E+02 7.2E-07 ];
INF_INVV                  (idx, [1:   4]) = [  5.49068E-08 0.00106  2.25961E-06 0.00045 ];

% Total scattering cross sections:

INF_SCATT0                (idx, [1:   4]) = [  5.34740E-01 0.00021  1.22311E+00 0.00054 ];
INF_SCATT1                (idx, [1:   4]) = [  2.36109E-01 0.00036  3.31622E-01 0.00109 ];
INF_SCATT2                (idx, [1:   4]) = [  9.36426E-02 0.00041  8.48904E-02 0.00239 ];
INF_SCATT3                (idx, [1:   4]) = [  7.42113E-03 0.00682  2.57589E-02 0.00763 ];
INF_SCATT4                (idx, [1:   4]) = [ -9.35660E-03 0.00512 -5.35691E-03 0.02804 ];
INF_SCATT5                (idx, [1:   4]) = [  5.50576E-04 0.07472  4.41052E-03 0.03628 ];
INF_SCATT6                (idx, [1:   4]) = [  5.08125E-03 0.00836 -1.21218E-02 0.01053 ];
INF_SCATT7                (idx, [1:   4]) = [  7.76475E-04 0.04598 -6.06858E-04 0.18504 ];

% Total scattering production cross sections:

INF_SCATTP0               (idx, [1:   4]) = [  5.34779E-01 0.00021  1.22311E+00 0.00054 ];
INF_SCATTP1               (idx, [1:   4]) = [  2.36111E-01 0.00037  3.31622E-01 0.00109 ];
INF_SCATTP2               (idx, [1:   4]) = [  9.36436E-02 0.00041  8.48904E-02 0.00239 ];
INF_SCATTP3               (idx, [1:   4]) = [  7.42065E-03 0.00681  2.57589E-02 0.00763 ];
INF_SCATTP4               (idx, [1:   4]) = [ -9.35670E-03 0.00511 -5.35691E-03 0.02804 ];
INF_SCATTP5               (idx, [1:   4]) = [  5.50394E-04 0.07490  4.41052E-03 0.03628 ];
INF_SCATTP6               (idx, [1:   4]) = [  5.08170E-03 0.00837 -1.21218E-02 0.01053 ];
INF_SCATTP7               (idx, [1:   4]) = [  7.76441E-04 0.04611 -6.06858E-04 0.18504 ];

% Diffusion parameters:

INF_TRANSPXS              (idx, [1:   4]) = [  2.30824E-01 0.00062  8.88636E-01 0.00065 ];
INF_DIFFCOEF              (idx, [1:   4]) = [  1.44411E+00 0.00062  3.75111E-01 0.00065 ];

% Reduced absoption and removal:

INF_RABSXS                (idx, [1:   4]) = [  1.13809E-02 0.00081  1.19401E-01 0.00092 ];
INF_REMXS                 (idx, [1:   4]) = [  2.72280E-02 0.00042  1.21819E-01 0.00116 ];

% Poison cross sections:

INF_I135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_YIELD          (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_I135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM147_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM148M_MICRO_ABS      (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_PM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_XE135_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_SM149_MACRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison decay constants:

PM147_LAMBDA              (idx, 1)        =  8.37272E-09 ;
PM148_LAMBDA              (idx, 1)        =  1.49451E-06 ;
PM148M_LAMBDA             (idx, 1)        =  1.95433E-07 ;
PM149_LAMBDA              (idx, 1)        =  3.62737E-06 ;
I135_LAMBDA               (idx, 1)        =  2.93061E-05 ;
XE135_LAMBDA              (idx, 1)        =  2.10657E-05 ;
XE135M_LAMBDA             (idx, 1)        =  7.55556E-04 ;
I135_BR                   (idx, 1)        =  8.34914E-01 ;

% Fission spectra:

INF_CHIT                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHIP                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
INF_CHID                  (idx, [1:   4]) = [  1.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

INF_S0                    (idx, [1:   8]) = [  5.18929E-01 0.00021  1.58107E-02 0.00079  2.33854E-03 0.00885  1.22077E+00 0.00054 ];
INF_S1                    (idx, [1:   8]) = [  2.31549E-01 0.00036  4.56063E-03 0.00202  9.24066E-04 0.01417  3.30698E-01 0.00107 ];
INF_S2                    (idx, [1:   8]) = [  9.50207E-02 0.00042 -1.37804E-03 0.00628  5.04758E-04 0.01895  8.43856E-02 0.00241 ];
INF_S3                    (idx, [1:   8]) = [  9.03606E-03 0.00520 -1.61493E-03 0.00529  1.66103E-04 0.04547  2.55928E-02 0.00765 ];
INF_S4                    (idx, [1:   8]) = [ -8.82410E-03 0.00521 -5.32499E-04 0.01108 -1.00244E-05 0.74167 -5.34689E-03 0.02778 ];
INF_S5                    (idx, [1:   8]) = [  5.27126E-04 0.07818  2.34504E-05 0.26788 -9.09433E-05 0.08001  4.50146E-03 0.03494 ];
INF_S6                    (idx, [1:   8]) = [  5.19935E-03 0.00822 -1.18097E-04 0.04014 -1.05804E-04 0.05580 -1.20160E-02 0.01036 ];
INF_S7                    (idx, [1:   8]) = [  9.19201E-04 0.03800 -1.42726E-04 0.03735 -9.84972E-05 0.07485 -5.08361E-04 0.22108 ];

% Scattering production matrixes:

INF_SP0                   (idx, [1:   8]) = [  5.18968E-01 0.00021  1.58107E-02 0.00079  2.33854E-03 0.00885  1.22077E+00 0.00054 ];
INF_SP1                   (idx, [1:   8]) = [  2.31550E-01 0.00036  4.56063E-03 0.00202  9.24066E-04 0.01417  3.30698E-01 0.00107 ];
INF_SP2                   (idx, [1:   8]) = [  9.50216E-02 0.00042 -1.37804E-03 0.00628  5.04758E-04 0.01895  8.43856E-02 0.00241 ];
INF_SP3                   (idx, [1:   8]) = [  9.03558E-03 0.00520 -1.61493E-03 0.00529  1.66103E-04 0.04547  2.55928E-02 0.00765 ];
INF_SP4                   (idx, [1:   8]) = [ -8.82420E-03 0.00520 -5.32499E-04 0.01108 -1.00244E-05 0.74167 -5.34689E-03 0.02778 ];
INF_SP5                   (idx, [1:   8]) = [  5.26944E-04 0.07837  2.34504E-05 0.26788 -9.09433E-05 0.08001  4.50146E-03 0.03494 ];
INF_SP6                   (idx, [1:   8]) = [  5.19980E-03 0.00823 -1.18097E-04 0.04014 -1.05804E-04 0.05580 -1.20160E-02 0.01036 ];
INF_SP7                   (idx, [1:   8]) = [  9.19167E-04 0.03811 -1.42726E-04 0.03735 -9.84972E-05 0.07485 -5.08361E-04 0.22108 ];

% Micro-group spectrum:

B1_MICRO_FLX              (idx, [1: 140]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Integral parameters:

B1_KINF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_KEFF                   (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_B2                     (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];
B1_ERR                    (idx, [1:   2]) = [  0.00000E+00 0.0E+00 ];

% Critical spectra in infinite geometry:

B1_FLX                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS_FLX               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reaction cross sections:

B1_TOT                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CAPT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_ABS                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_FISS                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NSF                    (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_NUBAR                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_KAPPA                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_INVV                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering cross sections:

B1_SCATT0                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT1                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT2                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT3                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT4                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT5                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT6                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATT7                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Total scattering production cross sections:

B1_SCATTP0                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP1                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP2                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP3                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP4                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP5                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP6                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SCATTP7                (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Diffusion parameters:

B1_TRANSPXS               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_DIFFCOEF               (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Reduced absoption and removal:

B1_RABSXS                 (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_REMXS                  (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Poison cross sections:

B1_I135_YIELD             (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_YIELD           (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_YIELD            (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_I135_MICRO_ABS         (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM147_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM148M_MICRO_ABS       (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_PM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MICRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_XE135_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SM149_MACRO_ABS        (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Fission spectra:

B1_CHIT                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHIP                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_CHID                   (idx, [1:   4]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering matrixes:

B1_S0                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S1                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S2                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S3                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S4                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S5                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S6                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_S7                     (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Scattering production matrixes:

B1_SP0                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP1                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP2                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP3                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP4                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP5                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP6                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];
B1_SP7                    (idx, [1:   8]) = [  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00  0.00000E+00 0.0E+00 ];

% Additional diffusion parameters:

CMM_TRANSPXS              (idx, [1:   4]) = [  2.42504E-01 0.00113  7.90689E-01 0.00592 ];
CMM_TRANSPXS_X            (idx, [1:   4]) = [  2.42313E-01 0.00172  7.95643E-01 0.00974 ];
CMM_TRANSPXS_Y            (idx, [1:   4]) = [  2.43122E-01 0.00216  7.97628E-01 0.01009 ];
CMM_TRANSPXS_Z            (idx, [1:   4]) = [  2.42131E-01 0.00229  7.82632E-01 0.01010 ];
CMM_DIFFCOEF              (idx, [1:   4]) = [  1.37459E+00 0.00113  4.21925E-01 0.00586 ];
CMM_DIFFCOEF_X            (idx, [1:   4]) = [  1.37573E+00 0.00172  4.19899E-01 0.00970 ];
CMM_DIFFCOEF_Y            (idx, [1:   4]) = [  1.37121E+00 0.00215  4.18919E-01 0.01001 ];
CMM_DIFFCOEF_Z            (idx, [1:   4]) = [  1.37684E+00 0.00230  4.26956E-01 0.01009 ];

% Delayed neutron parameters (Meulekamp method):

BETA_EFF                  (idx, [1:  18]) = [  7.07783E-03 0.01686  1.93251E-04 0.09911  1.05645E-03 0.04354  6.29916E-04 0.05654  1.33656E-03 0.03645  2.29965E-03 0.02864  6.94495E-04 0.05078  6.16753E-04 0.05472  2.50759E-04 0.07982 ];
LAMBDA                    (idx, [1:  18]) = [  4.63490E-01 0.02614  1.24667E-02 0.0E+00  2.82917E-02 0.0E+00  4.25244E-02 5.7E-09  1.33042E-01 0.0E+00  2.92467E-01 0.0E+00  6.66488E-01 0.0E+00  1.63478E+00 4.9E-09  3.55460E+00 0.0E+00 ];

