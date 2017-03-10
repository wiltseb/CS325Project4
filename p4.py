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
            if i != j:
                city2Identifyer = vertices[j][0]
                city2x = vertices[city2Identifyer][1]
                city2y = vertices[city2Identifyer][2]        
                G[city1Identifyer][city2Identifyer] = getDistance(city1x, city1y, city2x, city2y)
    return G

def sortEdges(adjacentVertices):
    '''
    sortedAdjList = []
    sortedGraph = {}
    for i in range(len(G)):
        currVertex = G[i]
        sortedAdjList = sorted(currVertex.items(), key=operator.itemgetter(1))
        sortedGraph[i] = sortedAdjList
    '''
    #https://www.ics.uci.edu/~eppstein/PADS/MinimumSpanningTree.py
    sortedGraph = sorted((G[u][v],u,v) for u in G for v in G[u])
    return sortedGraph
    



def primms_mst(G, startCityID):    
    for u in G:
        G[u]['distanceTo'] = float('inf')
    Gcopy = G.copy()
    currCityID = startCityID
    parent = {}
    while Gcopy:   
        u = Gcopy.pop(currCityID, None)
        currMin = float('inf')
        minVID = -1
        for v in u:
            if u[v] < currMin and v != 'distanceTo':
                currMin = u[v]
                minVID = v
                G[v]['distanceTo'] = u[v]
        parent[v] = currCityID
        currCityID = minVID
        print currCityID
        
                
            
    return parent

'''
FUNCTION USAGE IN END PRODUCT(Must accept problem instances on command line):
getInputData(sys.arv[1])
TSP function call
createOutputFile(tourLength, tour, sys.argv[2])
'''



vertices = getInputData("tsp_example_1.txt");
G = buildConnectedGraph(vertices)
primms_mst(G, 0)










    




