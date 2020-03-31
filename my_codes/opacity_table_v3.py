#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 11 15:25:45 2019

@author: ganguly93
"""

"""
This python script contains the class and
functions for the opacity table for the line
profiles project.
"""
import numpy as np
def bilinear_interp(q11, q21, q12, q22, x1, x2, y1, y2, x, y):
    x2x1 = x2 - x1
    y2y1 = y2 - y1
    x2x = x2 - x
    y2y = y2 - y
    yy1 = y - y1
    xx1 = x - x1
    return (q11*x2x*y2y + q21*xx1*y2y + q12*x2x*yy1 + q22*xx1*yy1)/(x2x1 * y2y1)

class opacity_table:
    def __init__(self,data_file):
        self.data = []
        """
        Read in the ion data.
        """
        print ("Reading the data file {}.".format(data_file))
        with open(data_file) as f:
            for line in f:
                self.data.append(line.split())

        print ("Done!")
        """
        break it down into smaller chunks
        """
        self.data = np.array(self.data[1:])
        self.xis = self.data[:,1].astype(float)
        self.temps = np.log10(self.data[:,0].astype(float)*1e4)
        self.data = np.array(self.data[:,2:])
        self.num_lines = len(self.data[0])/3
        self.num_temps = len(np.unique(self.temps))
        self.num_xis = len(np.unique(self.xis))
        
        line_arr = []
        
        for i in range(int(self.num_lines)):
            line_arr.append(self.data[0][1+i*3])

        self.lines = line_arr
        self.opact_dict = dict()

        print ("Creating interpolation functions for the following lines:")
        print (self.lines)

        for i,line in enumerate(line_arr):
            self.opact_dict.update({line :
                self.data[:,2+i*3].astype(float)})

    """
    Function to get the opacities after we've
    initialized the tables and such.
    we get opacities by passing log(xi) and log(T)
    """
    def get_opacity(self,xi,T,i,verbose=False):#,line=None):
        if i==0:
            line = self.lines[0]
        else:
            line = self.lines[1]
        l1 = np.where(np.unique(self.temps) < T)[0][-1]*self.num_xis
        l2 = np.where(np.unique(self.xis) < xi)[0][-1]
        l3 = l1 + self.num_xis
        q11 = self.opact_dict[line][l1+l2]
        q12 = self.opact_dict[line][l1+l2+1]
        q21 = self.opact_dict[line][l3+l2]
        q22 = self.opact_dict[line][l3+l2+1]
        x1 = self.temps[l1+l2]
        x2 = self.temps[l3+l2]
        y1 = self.xis[l1+l2]
        y2 = self.xis[l3+l2+1]
        ans = bilinear_interp(q11, q21, q12, q22, x1, x2, y1, y2, T, xi)

        if verbose:
            print (q11,q12,q21,q22,x1,x2,y1,y2)
            print (ans)
        return ans