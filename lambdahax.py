#!/bin/python2

import argparse
import csv
import os
import sys
import numpy as np

def isLast(itr):
    old = itr.next()
    for new in itr:
        yield False, old
        old = new
    yield True, old

os.system('clear')

maxEGO = 1.45
minEGO = 0.55
maxETE = 103
minCHT = 70
maxAbsDTPS = 0.0

parser = argparse.ArgumentParser(description='Hacks at FreeEMS CSV files')
parser.add_argument('logfile', nargs=1, help='FreeEMS CSV log file')

args = parser.parse_args()

if args.logfile == None:
    parser.print_help()
    sys.exit("\nNo logfile file. Exiting.")

logfilename = args.logfile[0]

tableAccWeight = np.genfromtxt("basetable.csv", delimiter=',')
tableAccLambda = np.genfromtxt("basetable.csv", delimiter=',')
tableLambda = np.genfromtxt("basetable.csv", delimiter=',')
np.set_printoptions(linewidth=160, precision=2, nanstr='', suppress=True)

with open(logfilename) as csvfile:
    dialect = csv.Sniffer().sniff(csvfile.read(1024))
    csvfile.seek(0)
    reader = csv.reader(csvfile, dialect, delimiter=',')

    row = reader.next()
    indexCHT = row.index("CHT")
    indexEGO = row.index("EGO")
    indexMAP = row.index("MAP")
    indexRPM = row.index("RPM")
    indexTPS = row.index("TPS")
    indexETE = row.index("ETE")

    row = reader.next()
    maxMAP = float(row[indexMAP])
    maxRPM = float(row[indexRPM])
    minMAP = float(row[indexMAP])
    minRPM = float(row[indexRPM])
    maxDTPS = 0

    lastTPS = 0
    datapoints = 0; ignoredEdge = 0; ignoredEGO = 0; ignoredETE = 0; ignoredDTPS=0;
    for (isLastDataPoint, row) in isLast(reader):
        datapoints += 1

        #tableAccWeight = np.genfromtxt("basetable.csv", delimiter=',')

        CHT  = float(row[indexCHT])
        EGO  = float(row[indexEGO])
        ETE  = float(row[indexETE])
        MAP  = float(row[indexMAP])
        RPM  = float(row[indexRPM])
        TPS  = float(row[indexTPS])
        DTPS = TPS - lastTPS
        lastTPS = TPS

        maxMAP = max(maxMAP, MAP)
        maxRPM = max(maxRPM, RPM)
        minMAP = min(minMAP, MAP)
        minRPM = min(minRPM, RPM)
        maxDTPS = max(maxDTPS, DTPS)

        r = 1; m = 1
        while r < tableAccWeight.shape[0]-1 and RPM > float(tableAccWeight[0,r]):
            r += 1
        while m < tableAccWeight.shape[1]-1 and MAP > float(tableAccWeight[m,0]):
            m += 1

        if ETE > maxETE: # skip warm-up
            ignoredETE += 1
        elif abs(DTPS) > maxAbsDTPS: # skip accel transients
            ignoredDTPS += 1
        elif EGO >= maxEGO or EGO <= minEGO: # skip saturated EGO readings
            ignoredEGO += 1
        elif r >= 11 or m >= 11 or r <= 1 or m <= 1: # TODO fix the edge cases - just skip for now
            ignoredEdge += 1
        else:
            mapH = float(tableAccWeight[m,0])
            mapL = float(tableAccWeight[m-1,0])
            rpmH = float(tableAccWeight[0,r])
            rpmL = float(tableAccWeight[0,r-1])

            assert mapH >= mapL, "%f >= %f" % (mapH, mapL)
            assert rpmH >= rpmL, "%f >= %f" % (rpmH, rpmL)
            assert MAP >= mapL,  "%f >= %f" % (MAP, mapL)
            assert RPM >= rpmL,  "%f >= %f" % (RPM, rpmL)
            assert mapH >= MAP,  "%f >= %f" % (mapH, MAP)
            assert rpmH >= RPM,  "%f >= %f" % (rpmH, RPM)

            cell = (rpmH-rpmL) * (mapH-mapL)
            assert cell > 0

            m1r1 = (rpmH-RPM)*(mapH-MAP) / cell
            m1r0 = (RPM-rpmL)*(mapH-MAP) / cell
            m0r1 = (rpmH-RPM)*(MAP-mapL) / cell
            m0r0 = (RPM-rpmL)*(MAP-mapL) / cell
            assert abs(1 - m1r1 - m1r0 - m0r1 - m0r0) < 0.0001

            tableAccWeight[m-1,r-1] += m1r1
            tableAccWeight[m-1,r]   += m1r0
            tableAccWeight[m,r-1]   += m0r1
            tableAccWeight[m,r]     += m0r0

            tableAccLambda[m-1,r-1] += m1r1 * EGO
            tableAccLambda[m-1,r]   += m1r0 * EGO
            tableAccLambda[m,r-1]   += m0r1 * EGO
            tableAccLambda[m,r]     += m0r0 * EGO

            tableLambda[m-1,r-1] = tableAccLambda[m-1,r-1] / tableAccWeight[m-1,r-1]
            tableLambda[m-1,r]   = tableAccLambda[m-1,r]   / tableAccWeight[m-1,r]
            tableLambda[m,r-1]   = tableAccLambda[m,r-1]   / tableAccWeight[m,r-1]
            tableLambda[m,r]     = tableAccLambda[m,r]     / tableAccWeight[m,r]

        if (datapoints % 1000 == 0) or isLastDataPoint:
            ignoredPoints = ignoredEdge + ignoredEGO + ignoredETE + ignoredDTPS
            if isLastDataPoint:
                os.system('clear')

            n = 80; sys.stdout.write('\033[{}A'.format(n))
            print 'Cell Weight'
            print tableAccWeight
            print '\nAccumulated EGO'
            print tableAccLambda
            print '\nMean Lambda'
            print tableLambda
            print '\nMAP: {:3.1f}'.format(MAP), '\tRPM: {:4.0f}'.format(RPM), '\tEGO: {:1.2f}'.format(EGO), '\tTPS: {:1.2f}'.format(TPS), '\tDTPS: {:1.2f}'.format(DTPS), '\tCHT: {:1.2f}'.format(CHT), '\tETE: {:1.2f}'.format(ETE)
            print
            print 'minRPM:', minRPM, '     '
            print 'maxRPM:', maxRPM, '     '
            print 'minMAP:', minMAP, '     '
            print 'maxMAP:', maxMAP, '     '
            print 'maxDTPS:', maxDTPS, '     '
            print
            print 'Ignored', ignoredPoints, 'of', datapoints, 'data points -', 100*ignoredPoints/datapoints, '% waste   '
            print '\t- ETE =', ignoredETE, 100*ignoredETE/ignoredPoints, '%   '
            print '\t- DTPS =', ignoredDTPS, 100*ignoredDTPS/ignoredPoints, '%   '
            print '\t- EGO =', ignoredEGO, 100*ignoredEGO/ignoredPoints, '%   '
            print '\t- Edge =', ignoredEdge, 100*ignoredEdge/ignoredPoints, '%   '



