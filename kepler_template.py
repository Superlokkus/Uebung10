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
    return opt.newton(KeplerEq,3, args=(Planet[planet]["numEx"],t,t0,Planet[planet]["Periode"]))

# Weiter mit Berechnung der wahren Anomalie und des Abstandes von der Sonne:

def W_Anom(planet,t=0,t0=0):
    """Wahre Anomalie"""
    E = ex_Anom(planet,t,t0)
    eps = Planet[planet]["numEx"] #Zur besseren Lesbarkeit
    
    if (E >= 0 and E <= sc.pi):
        return sc.arccos((sc.cos(E) - eps) / (1 - eps*sc.cos(E)))
        
    elif (E >= sc.pi and E <= 2*sc.pi):
        return 2*sc.pi - sc.arccos((sc.cos(E) - eps) / (1 - eps*sc.cos(E)))
            
    else: 
        return None #Außerhalb der Intervalle nicht definiert
        
def AbstandSonne(planet,t=0,t0=0):
    a = Planet[planet]["Halbachse"] #Zur besseren Lesbarkeit
    eps = Planet[planet]["numEx"] #Zur besseren Lesbarkeit
    
    return a * (1 - eps**2) / (1 + eps*sc.cos(W_Anom(planet,t,t0)))

print W_Anom("Erde",300,0)
#print sc.linspace(0,Planet["Erde"]["Periode"],1000)
#print W_Anom("Erde",sc.linspace(0,Planet["Erde"]["Periode"],1000))

plt.polar(W_Anom("Erde",sc.linspace(0,Planet["Erde"]["Periode"],1000)),AbstandSonne("Erde"),sc.linspace(0,Planet["Erde"]["Periode"],1000))

plt.show()
