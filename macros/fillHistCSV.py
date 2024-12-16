#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *
#import cmsstyle as CMS

import subprocess
import sys
import os
import random

# uses python3
# jwk code


def rolldie6() :

    die6 = [1,2,3,4,5,6]
    return random.choice(die6)

def explode( rolls, ndie ) :
    
    #print( rolls )
    if ndie == 1 :
        return False
    for num in rolls :
        if num != 6 : 
            return False
    return True

def rollndie( num, rolls ) :

    count = 0
    while count < num :
        rolls.append( rolldie6() )
        count = count + 1
    
def makeroll( ndie ) :

    rolls = []
    rollndie( ndie, rolls )
    excnt = 0;
    #print( '---------------------')
    while explode( rolls, ndie ) :
        rollndie( ndie, rolls )
        excnt = excnt + 1
        #if excnt > 1 : 
        #    print( "explode ", excnt )
        if excnt > 100 : break
    rolls.sort()
    #print( rolls )
    rcnt = 0
    rsum = 0
    fndie = len(rolls)
    for roll in rolls :
        rsum = rsum + roll
        if rsum < 10 : rcnt = rcnt + 1
    if rsum == ndie : return 0
    if rsum < 10 : return 0
    return fndie - rcnt

def makeHist() :

    outfilename = "csvfile.root"
    mtfile = TFile.Open(outfilename,"RECREATE")
    mtfile.cd()
    h2d = TH2D("h2d","probplot;nDice;nSucesses",21,-0.5,20.5,21,-0.5,20.5)
    res = [0] * 20
    dis = []
    #print( res )
    for ndie in range(1,20) :
        for trail in range(1,1000000) :
            wins = makeroll( ndie )
            #print( "wins: ", wins )
            if wins > 20 : 
                wins = 20
            res[wins] += 1
        #print( res )
        for suc in range( 0, 20 ) :
            dis.append(float(res[suc])/999999.0)
            h2d.Fill( ndie, suc, float(res[suc])/999999.0 )
            res[suc] = 0
        print(dis)
        dis.clear()
    h2d.Write()
    

def myReadFile() :

    outfilename = "csvfile.root"
    mtfile = TFile.Open(outfilename,"RECREATE")
    mtfile.cd()

    col = 0
    h2d = TH2D("h2d","probplot",12,0,12,20,0,20)
    myinfile = open("grid8.txt","r")
    thelines = myinfile.readlines()   
    #print( thelines )
 
    for theline in thelines :

        print( theline )
        row = 0
        data = theline.split(",")
        for entry in data :
            print( entry )
            if entry : prob = float( entry ) 
            else : prob = 0
            #print( "col: ", col, "Row: ", row, " data: ",prob)
            h2d.Fill( col, row, prob )
            row = row + 1
        col = col + 1
       
    h2d.Write()
    mtfile.Close()

makeHist()
#myReadFile()
