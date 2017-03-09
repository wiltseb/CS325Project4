import math
import sys

'''
Takes in input from a file with the following format:
Each line defines a city and each line has three numbers separated by
white space.
o The first number is the city identifier
o The second number is the city's x-coordinate
o The third number is the city's y-coordinate.

Returns a list of lists in the format:
[identifyer, x-coordinate, y-coordinate]
'''
def getInputData(filename):
    data = []
    with open(filename) as inFile:
        for line in inFile:
            line = line.split()
            line = [int(i) for i in line]
            data.append(line)
    return data


'''
Gets distance between ordered pairs (x1,y1) and (x2,y2) rounded to the
nearest integer.
'''
def getLocalDistance(x1,y1,x2,y2):
    x = pow((x1-x2), 2)
    y = pow((y1-y2), 2)
    return int(math.sqrt(x+y))

'''
Creates output file formatted to specification

params --
tourLength: total length of the tour (integer)
tour: list of lists in the form [city id, x-coordinate, y-coordinate] in order of the tour
outFilename: name of output file
'''
def createOutputFile(tourLength, tour, outFilename):
    outFile = open(outFilename, 'w')
    outFile.write(str(tourLength) + '\n')
    for city in tour:
        outFile.write(str(city[0]) + '\n')

'''
quick test to ensure that createOutputFile formats correctly
'''
def testCreateOutputFile():
    lst = [[0, 200, 800],[1, 3600, 2300],[2, 3100, 3300],[3, 4700, 5750]]    
    createOutputFile(500, lst, "testout.txt")


'''
FUNCTION USAGE IN END PRODUCT(Must accept problem instances on command line):
getInputData(sys.arv[1])
TSP function call
createOutputFile(tourLength, tour, sys.argv[2])
'''




