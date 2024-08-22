# Python3 script to call and test the NXT class
# Author: R. Guasch
# Date : 24 July 2024
# Dependencies : GEOM.py, CARCEL.py, NXT_CLUSTER.py
#
from GEOM import GEO
from CARCEL import CARCEL
from CAR2D import CAR2D
from NXT_CLUSTER import NXTGAC

import numpy as np


# Create a 2D cartesian cell geometry, at the first level
                                                              # 0.0, Rcomb1,   Rcomb2,   Rcomb3,  Rcomb4, Rgap,   Rclad
AT10_C1 = CARCEL("C1", level = 1, nr = 6, radii = np.array([0.0, 0.313602, 0.396678, 0.43227, 0.4435, 0.4520, 0.5140]), meshx = np.array([0.0, 1.295]), meshy = np.array([0.0, 1.295]))


# Analyse it with NXT_CLUSTER
NXT_C1 = NXTGAC(AT10_C1)