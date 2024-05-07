# Python3 script to test functionnalities of PyNE package 
# These involve the MCNP API and ACE file parser API
# Author : R. Guasch
# Part of the PyNjoy2016 XS survey
# Inteded use is to test the meta-stable isotopes ace files with PyNE parsing

import pyne

lib_Am242m_file = pyne.ace.Library('home/p117902/Serpent2/xs/jeff311/acedata_pynjoy2016/95342_550K.ace')

lib_Am242m_file.read()