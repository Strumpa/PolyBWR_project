#
# test2_lcm: non regression testing for lcm class with (de)serialization
#
import lcm
from assertS import *
import numpy as np
my_lcm=lcm.new('LCM','FILE1',impx=0)
ia=np.array([8, 7, 8, 4, 9, 1, 0, 4], dtype='i')
ra=np.array([8.0,6.0,5.0,2.0,1.0], dtype='f')
da=np.array([8.0,6.0,5.0,2.0,1.0], dtype='d')
my_lcm['key1']='new comments for this record'
my_lcm['key2']=ia
my_lcm['key3']=ra
my_lcm['key4']=da

# Export (serialize) the initial LCM file to ascii
print("****** Export the initial LCM file to ascii ******")
my_lcm2 = lcm.new('ASCII', 'EXPORT1', pyobj=my_lcm, impx=0)

# Reimport (deserialize) the LCM file
my_lcm3 = lcm.new('LCM_INP', 'EXPORT1', impx=0)
# Modify a record
ia=np.array([1, 1, 2, 2, 3, 3, 4, 4], dtype='i')
my_lcm3['key2']=ia

# Export (serialize) the modified LCM file to ascii
print("****** Export the modified LCM file to ascii ******")
my_lcm4 = lcm.new('ASCII', 'EXPORT2', pyobj=my_lcm3, impx=0)

print("test test2_lcm completed")
