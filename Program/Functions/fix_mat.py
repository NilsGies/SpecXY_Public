#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Nils B. Gies
"""
import scipy.io as sio

def fix_mat(f):
    struct=sio.loadmat(str(f))
    sio.savemat(str(f),struct)

fix_mat(input)