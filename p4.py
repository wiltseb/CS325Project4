import math
import sys
import operator


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
def getDistance(x1,y1,x2,y2):

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
Takes a list of vertices and builds a complete graph in the form:

G = {
    {cityID0:
         adjacentID1: dist(city0, city1),
         adjacentID2: dist(city0, city2),
        ...
         adjacentIDn: dist(city0, cityn)}
    {cityID1:
        {adjacentID0: dist(city1, city0),
         adjacentID2: dist(city1, city2),
        ...
         adjacentIDn: dist(city1, cityn)}
        }
        ....

So, G[cityA][cityB} returns the distance between city A and city B

'''
def buildConnectedGraph(vertices): 
    G = {}
    #put each vertex in adjacency list for every other vertex
    for i in range(len(vertices)):
        city1Identifyer = vertices[i][0]
        city1x = vertices[city1Identifyer][1]
        city1y = vertices[city1Identifyer][2]
        G[city1Identifyer] = {}
        for j in range(len(vertices)):
            if i != j:  #don't put city in its own adjacency list
                city2Identifyer = vertices[j][0]
                city2x = vertices[city2Identifyer][1]
                city2y = vertices[city2Identifyer][2]
                #G[u][v] is the distance between u and v
                G[city1Identifyer][city2Identifyer] = getDistance(city1x, city1y, city2x, city2y)
    return G


'''
Takes in a graph, G (a dictionary of dictionaries), and returns a new graph with keys "predecessor" and "distanceTo",
which indicate the current vertex's parent vertex and how far the current vertex and parent are away from each other, respectively
'''
def prims_mst(G, startCityID):
    #initialize each vertex in G with an infinite distanceTo
    keys = G.keys()
    for u in G:
        G[u]['distanceTo'] = float('inf')
    newG = {}
    currCityID = startCityID
    #Start city's predecessor is None, the distanceTo is 0
    G[currCityID]['predecessor'] = None;
    G[currCityID]['distanceTo'] = 0;
    #Take the city(and its adjacency list) with the closest distance out of original graph and put it in a new
    #graph with predecessor and distanceTo data
    while G:
        u = G.pop(currCityID, None)
        newG[currCityID] = u
        currMin = float('inf')
        for v in u:
            #Only interested in vertices still in G ('popped' vertices have already been visited)
            if v in G and v != 'predecessor' and v != 'distanceTo':
                if u[v] < G[v]['distanceTo'] :
                    G[v]['distanceTo'] = u[v]
                    G[v]['predecessor'] = currCityID
        minDist = float('inf')
        for i in keys:
            if i != 'predecessor' and i != 'distanceTo' and i in G:
                if G[i]['distanceTo'] < minDist:
                    minDist = G[i]['distanceTo']
                    currCityID = i
        
    return newG

'''
Takes in a graph that has MST data and prints MST heirarchy
Also prints number of cities for a sanity check
'''
def printMSTHeirarchy(MSTGraph):
    cityCount = 0
    for i in MSTGraph:
        print str(i) + "'s predecessor is: " + str(MSTGraph[i]['predecessor']) + ". Distance is: " + str(MSTGraph[i]['distanceTo'])
        cityCount += 1
    print "Number of cities: " + str(cityCount)
    totalEdgeWeights=0
    for i in MSTGraph:
        totalEdgeWeights += MSTGraph[i]['distanceTo']
    print "Total edge weights: " + str(totalEdgeWeights)
    tree = getMSTree(MSTGraph)
    print "Current Tree Structure: "
    print tree

'''
Takes in a graph with MST data (predecessor and distanceTo) and returns a tree only MST edges
'''
def getMSTree(MSTGraph):
    cityCount=0
    MSTree = {}
    for i in MSTGraph:
        adjCityID = MSTGraph[i]['predecessor']
        adjCityDist = MSTGraph[i]['distanceTo']
        if adjCityID != None:
            if i not in MSTree:
                MSTree[i] = {adjCityID: adjCityDist}
            else:
                MSTree[i][adjCityID] = adjCityDist
            if adjCityID not in MSTree:
                MSTree[adjCityID] = {i: adjCityDist}
            else:
                MSTree[adjCityID][i] = adjCityDist
    return MSTree

'''
Verified Prim's works with the visualization from class at:
http://rosulek.github.io/vamonos/demos/prims.html
'''
def testPrims():
    G = { 'a':{'b':6,'c':3,'e':9},'b':{'a':6,'c':4,'d':2,'h':5,'i':4},'c':{'a':3,'b':4,'d':2,'e':9},
         'd':{'b':2,'c':2,'f':8},'e':{'a':9,'c':9,'f':8,'g':18},'f':{'e':8,'d':8,'h':9,'g':10},
         'g':{'i':4,'h':3,'f':10,'e':18},'h':{'f':9,'g':3,'b':5},'i':{'b':4,'g':4}}
    printMSTHeirarchy(prims_mst(G, 'b'))



'''
FUNCTION USAGE IN END PRODUCT(Must accept problem instances on command line):
getInputData(sys.arv[1])
TSP function call
createOutputFile(tourLength, tour, sys.argv[2])
'''


'''
vertices = getInputData("tsp_example_1.txt");
G = buildConnectedGraph(vertices)
newG = prims_mst(G, 75)
'''

testPrims()










    




