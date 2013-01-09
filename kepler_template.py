#!/usr/bin/python
# encoding=utf-8

"""Übung 10 Numerische Nullstellensuche: Keplerproblem"""
#Markus Klemm WS12/13 Phy-BA


from __future__ import division  # problemlose Ganzzahl-Division
import scipy as sc               # scipy beinhaltet numpy
import scipy.optimize as opt     # Funktionen für Nullstellensuche
from matplotlib import pyplot as plt

"""
Nullstellensuche: Kepler-Problem
"""
erdjahr = 365 # Tage
# Parameter der Himmelskörper (Periode [d], mumerische Exzentrizität und
# große Halbachse [AE]):
Planet = {"Erde": {"Periode": erdjahr, "numEx": 0.0167, "Halbachse": 1.0}, 
          "Pluto": {"Periode": 248*erdjahr, "numEx": 0.25, "Halbachse": 39.44},  
          "Halley": {"Periode": 76*erdjahr, "numEx": 0.97, "Halbachse":17.94} }

def KeplerEq(E, eps, t, t0, T):
    """ Kepler-Gleichung """
    return E - eps*sc.sin(E) - 2*sc.pi*(t - t0)/T

def KeplerEq_prime(E, eps, t, t0, T):
    """ 1. Ableitung der Kepler-Gleichung """
    return 1 - eps*sc.cos(E)

def ex_Anom(planet, t=0, t0=0):
    """ Exzentrische Anomalie als Nullstelle der Kepler-Gleichung """
    return opt.newton(KeplerEq,3, args=(Planet["Erde"]["numEx"],t,t0,Planet["Erde"]["Periode"]))

# Weiter mit Berechnung der wahren Anomalie und des Abstandes von der Sonne:

def W_Anom(E,eps):
    """Wahre Anomalie"""
    if (E >= 0 and E <= sc.pi):
        return sc.arccos*((sc.cos(E) - eps) / (1 - eps*sc.cos(E)))
        
    elif (E >= sc.pi and E <= 2*sc.pi):
        return 2*sc.pi - sc.arccos*((sc.cos(E) - eps) / (1 - eps*sc.cos(E)))
            
    else: 
        return None
        
def AbstandSonne(t,a,eps):
    return a * (1 - eps**2) / (1 + eps*sc.cos(W_Anom()))


print ex_Anom(["Erde"])
