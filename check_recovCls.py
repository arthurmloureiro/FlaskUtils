#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Plots the RecovCls and the relative error betweem theory and signal C(l)s

Arthur Loureiro (UCL) - Jun/2019
"""
from __future__ import division, print_function
import sys
import numpy as np
import pylab as plt
import matplotlib as mpl
import os
import glob
from scipy.interpolate import InterpolatedUnivariateSpline as intp
from matplotlib import gridspec


def create_output_folder(config):
    """
    Creates an output folder for the plots.
    If folder already exists, creates a new one.
    
    Returns the name of the output folder!
    """
    pathtoConfig = os.path.dirname(config)
    rootName = os.path.splitext(os.path.basename(config))[0]
    
    newFolderName = "/check_RecovCls-"+rootName
    newFolder = os.path.join(pathtoConfig, newFolderName)
    
    if os.path.isdir(newFolder) == True:
        print("Folder ", newFolder, " exists.")
    else:
        try:
            print("Creating folder ", newFolder, " for the outputs.")
            os.mkdir(newFolder)
            
        except:
            print("Cannot create folder at ", pathtoConfig, "! \nWill create folder in current directory!")
            newFolder = "."+newFolderName
            if os.path.isdir(newFolder) == True:
                print("Folder ", newFolder, " exists.")
            else:
                os.mkdir(newFolder)
    
    return newFolder  
    

def search_flask_args(pathFlaskFile, argName):
    """
    returns the argument used in the config file for flask
    """
    with open(pathFlaskFile) as search:
        for line in search:
            #line = line.rstrip().replace(" ", "")
            if argName in line: #FIXME: need to find exact match!
                args = line.split(':')[1].rstrip().replace(" ", "").replace("\t", "").split('#')[0]
                break # since I can't find a way to find exact match
    
    if args[0]=='0':
        print("Argument not given to flask!")
        return None
    else:
        return args

def get_config_path(config):
    """
    Gets the flask config path
    """
    dirname = os.path.dirname(config)
    return dirname

def get_recov_cls_dict(config):
    """
    Returns the RecovCls in an array
    """
    
    configDirName = get_config_path(config)
    
    recovClsPathFromConfig = search_flask_args(config, "RECOVCLS_OUT")
    
    if recovClsPathFromConfig is None:
        print("No recovCls were calculated by Flask. \nDiagnosis cannot be performed!\nExiting the script...")
        sys.exit(-1)
    else:
        #check if the path is absolute:
        if recovClsPathFromConfig[0] == '/':
            RecovClsPath = recovClsPathFromConfig
        else:
            RecovClsPath = os.path.join(configDirName, recovClsPathFromConfig)
        
        assert os.path.isfile(RecovClsPath), "Cannot find the config file!"
        
        RecovCls = np.loadtxt(RecovClsPath)
        
        # getting the header information to create the dictionary:
        with open(RecovClsPath) as header:
            for line in header:
                names = line.replace("Cl-", "").rstrip().split(' ')
                names.remove('#')
                break
        #RecovCls Dictionary
        DictRecovCls = {}

        for n, i in zip(names, RecovCls.T):
            DictRecovCls[n] = i
        
    return DictRecovCls


def get_input_cls_dict(config, RecovDict):
    """
    recovers just the list of files!
    """
    configDirName = get_config_path(config)
    
    ClsPath = search_flask_args(config, "CL_PREFIX")
    
    if ClsPath is None:
        print("No input cls in the Flask config! \nDiagnosis cannot be performed!\nExiting the script...")
        sys.exit(-1)
    else:
        #check if the path is absolute:
        if ClsPath[0] == '/':
            ClsPath = ClsPath
        else:
            ClsPath = os.path.join(configDirName, ClsPath)
            
        InputClsDict = {}
        for k in RecovDict.keys():
            FileCl = glob.glob(ClsPath + k + ".dat")
            
            if FileCl != []:
                print("Reading input cls from: ", FileCl)
                assert os.path.isfile(FileCl[0]), "Cannot find the file for input cls"
                ell, cls = np.loadtxt(FileCl[0], unpack=1)
                InputClsDict['l'] = ell
                InputClsDict[k] = cls
        
        return InputClsDict

def check_ell_range(recovcls, inputcls):
    """
    check if ell ranges are the same for input and recov:
    True = Same Ell range
    False = Different Ell range
    """
    
    if (recovcls['l'].min() == inputcls['l'].min()) and (recovcls['l'].max() == inputcls['l'].max()):
        return True
    else:
        print("Different ell ranges.") 
        print("\nRecov: ", recovcls['l'].min(), " ", recovcls['l'].max())
        print("\nInput: ", inputcls['l'].min(), " ", inputcls['l'].max())
        return False

def get_nside(config):
    """
    returns the nside
    """
    return int(search_flask_args(config, "NSIDE"))
        
        
    
    
def plot_recov_vs_input(recovcls, inputcls, outputFolder):
    """
    Plots all of the input cls compared with the recov ones
    Plots also the relative % error between them
    
    FIXME: PIXEL WINDOW FUNCTION!
    """
    
    if check_ell_range(recovcls,inputcls) == False:
        ell_max = min(recovcls['l'].max(), inputcls['l'].max())
        ell_min = max(recovcls['l'].min(), inputcls['l'].min())
        print("Using the following range: ", ell_min, ell_max)
    else:
        ell_min = recovcls['l'].min()
        ell_max = inputcls['l'].max()
    
    ell_vec = np.arange(ell_min, ell_max + 1)
    fact = ell_vec*(ell_vec + 1)
    
    # checking if pixwin function was applied:
    if search_flask_args(config, "APPLY_PIXWIN") == '1':
        import healpy as hp
        print("Flask simulation using pixel window function. Will apply corrections...")
        nside = get_nside(config)
        pixWin2 = (hp.pixwin(nside)[int(ell_min):(int(ell_max) + 1)])**2
    else:
        pixWin2 = np.ones_like(ell_vec)
    
    for k in recovcls.keys():
        #interpolates both functions because that's what we have:
        if k != 'l':
            splRecov = intp(recovcls['l'], recovcls[k])
            splInput = intp(inputcls['l'], inputcls[k])
            
            fig = plt.figure()
            gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])
            plt.subplot(gs[0])
            plt.title(k + "\n(f1 = Pos; f2 = shear)")
            plt.loglog()
            plt.ylabel(r"$\ell(\ell +1)C_{\ell}$")
            plt.xlabel(r"$\ell$")
            plt.plot(ell_vec, fact*splRecov(ell_vec)/pixWin2, label="Recov")
            plt.plot(ell_vec, fact*splInput(ell_vec), label="Input")
            plt.legend(loc=0)
            
            plt.subplot(gs[1])
            plt.ylabel(r"Frac. Error $\%$")
            plt.xlabel(r"$\ell$")
            #plt.xscale('log')
            plt.ylim(-15,15)
            plt.plot(ell_vec, (splInput(ell_vec) / (splRecov(ell_vec)/pixWin2) - 1)*100, label="frac error")
            plt.axhline(0, ls='--', lw=2.0)
            
            
            plt.show()
    
    return None

def main(config):
    """
    Main part of the script
    """
    
    # uncomment the following line if you have no display access
    #plt.switch_backend('agg')
    
    assert os.path.isfile(config), "Cannot find the config file!"
    
    # tries to create a folder in the same dir as the .config
    # if it fails, creates in the current directory
    outputFolder = create_output_folder(config) 
    
    recovClsDict = get_recov_cls_dict(config)
    
    inputClsDict =  get_input_cls_dict(config, recovClsDict) #FIXME: RECOV CLS FORMAT DON'T MACH THIS LIST!
    
    plot_recov_vs_input(recovClsDict, inputClsDict, outputFolder)
    # print(recovClsDict.keys())
    # print("\n")
    # print(inputClsDict.keys())
    
    
if __name__=='__main__':
    if sys.argv[1:] == ["-h"]:
        print("USAGE: ./check_recovCls.py </path/to/flask.config> \nPlots Flask RecovCls compared to the input theory. \nPlots can be found in a new folder called 'check_RecovCls-<configname>' in the same directory as the .config.")
        sys.exit(-1)
    else:
        config = sys.argv[1]
    main(config)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    