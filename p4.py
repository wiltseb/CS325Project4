import math
import sys
import operator
import copy
from collections import defaultdict
import time


MAX_TIME = 180

'''
Graph for testing purposes, adapted from Prim's Visualization class module:
http://rosulek.github.io/vamonos/demos/prims.html

Note this one is not complete like the graphs
for TSP will be, but for functions for MST and reducing G it should work the same way.
'''

testGraph = { 'a':{'b':6,'c':3,'e':9},'b':{'a':6,'c':4,'d':2,'h':5,'i':4},'c':{'a':3,'b':4,'d':2,'e':9},
 'd':{'b':2,'c':2,'f':8},'e':{'a':9,'c':9,'f':8,'g':18},'f':{'e':8,'d':8,'h':9,'g':10},
 'g':{'i':4,'h':3,'f':10,'e':18},'h':{'f':9,'g':3,'b':5},'i':{'b':4,'g':4}}


'''
NEAREST NEIGHBOR FUNCTIONS
based on pseudocode from
https://en.wikipedia.org/wiki/Nearest_neighbour_algorithm
'''


'''
Takes in a city Identifyer (u) and a list of cities in the form
[ [int(cityID), int(xPos), int(yPos), bool(visited)], [...], ... ]
returns the nearest city in the form:
[cityID, distanceFrom(u)]
'''
def getNearest(u, cities):
    minDist = float('inf')
    nearestNeighbor = -1

    for v in range(len(cities)):
        if u != v and cities[v][3] == False:
            dist = getDistanceHelper(cities[u], cities[v])
            if dist < minDist:
                nearestNeighbor = v
                minDist = dist

    cities[nearestNeighbor][3] = True
    return [nearestNeighbor, minDist]

'''
Takes in a list of cities
Returns a list with the first index the TSP tour;
the next index is the total distance of the tour.
'''
def nearestNeighbor(cities, originID):
    TSPTour = [] #final tour

    #get first nearest neighbor
    first = getNearest(originID, cities)
    TSPTour.append(first[0])
    totalDist = first[1]

    nextCity = first
    #continue to find nearest neighbor iteratively
    for i in range(len(cities)-1):
        currCity = getNearest(nextCity[0], cities)
        TSPTour.append(currCity[0])
        totalDist += currCity[1]
        currCity = nextCity

    #last city in tour to origin
    totalDist += getDistanceHelper(cities[TSPTour[0]], cities[TSPTour[-1]])

    return [TSPTour, totalDist]





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
            line.append(False)
            data.append(line)
    return data

def getDistanceHelper(city1, city2):
    return getDistance(city1[1], city1[2], city2[1], city2[2])

'''
Gets distance between ordered pairs (x1,y1) and (x2,y2) rounded to the
nearest integer.
'''
def getDistance(x1,y1,x2,y2):
    dx = x2-x1
    dy = y2 - y1
    
    return int(round(math.sqrt(dx*dx + dy*dy)))

'''
Creates output file formatted to specification

params --
tourLength: total length of the tour (integer)
tour: list city IDs
outFilename: name of output file
'''
def createOutputFile(tour, tourLength, inFilename):
    outFilename = inFilename + ".tour"
    outFile = open(outFilename, 'w')
    outFile.write(str(tourLength) + '\n')
    for i in range(len(tour)-1):
        outFile.write(str(tour[i]) + '\n')
    outFile.write(str(tour[-1]) + '\n')

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
def buildCompleteGraph(vertices):
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
This Euler Tour is based on the algorithm pseudocode given here:
http://www.algorithmist.com/index.php/Euler_tour
The algorithm will find the Euler Tour recursively
'''


def eulerTour(multiGraph, tour, city):
    for neighbor in multiGraph[city]:
        if multiGraph[city][neighbor] > 0:
            multiGraph[city][neighbor] -= 1
            multiGraph[neighbor][city] -= 1
            eulerTour(multiGraph, tour, neighbor)
    tour.append(city)



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

'''
def testMSTReduce(Graph, startCityID):
    initialNumVerts = len(Graph)  #keep track of how many vertices we start with
    print "Original Graph:"
    for i in Graph:
        print "City " + str(i) + " adjacency list: " + str(Graph[i])
    print

    #return a graph that has data for a MST in G. All of newG's data is the same (same adjacency list)
    #as G, except now each city has a 'predecessor' and 'distanceTo'
    #Note, passing a copy of G because prims_mst modifies the graph passed in
    GraphCopy = copy.deepcopy(Graph)
    newG = prims_mst(GraphCopy, startCityID)
    assert( len(newG) == initialNumVerts ) #Shouldn't lose any vertices

    print "Original Graph with MST data:"
    for i in newG:
        print "City " + str(i) + " adjacency list: " + str(newG[i])
    print

    #Return a tree from G such that the only connections are those that make up the MST
    tree = getMSTree(newG)
    assert( len(tree) == initialNumVerts ) # shouldn't lose any vertices
    print "The MST in G (only connections are those from MST:"
    for i in tree:
        print "City " + str(i) + " adjacency list: " + str(tree[i])
    print


    #Return a list of cityIdentifyers (vertices) that have odd degree
    vwod = getVerticesWithOddDegree(tree)
    numVWOD = len(vwod)  #keep track of how many vertices with odd degree
    print "Vertices in MStree with odd degree:"
    print vwod
    print

    #Reduce the original graph, G, to have only vertices of odd degree
    GPrime = reduceG(Graph, vwod)
    assert( len(GPrime) == numVWOD ) #GPrime should have only vertices with odd degree
    print "Original graph with all original connections, reduced to have only MST vertices with odd degree:"
    for i in GPrime:
        print "City " + str(i) + " adjacency list: " + str(GPrime[i])
    print

    #Find perfect matching with minimum weight
    matching = greedyMatching(GPrime, vwod)
    print "Perfect matching with minimum weight from MST vertices with odd degree:"
    for i in matching:
        cityOne = "From city {}".format(i)
        cityTwo = "to city {}".format(matching[i].keys()[0])
        weight = "with weight {}".format(matching[i].values()[0])
        print "{} {} {}".format(cityOne, cityTwo, weight)


'''
Input: MStree
Output: list of vertices with odd degree from MST
'''
def getVerticesWithOddDegree(MSTree):
    vwod = []
    for i in MSTree:
        if len(MSTree[i]) % 2 == 1:
            vwod.append(i)
    return vwod

'''
Reduce original graph to have only vertices from MST with odd degree
'''
def reduceG(originalGraph, vwod):
    reducedG = {}
    for i in vwod:
        reducedG[i] = originalGraph[i]
    return reducedG



def greedyMatching(oddGraph, vwod):
    '''
    Greedy algorithm that finds the lowest-weight matching edges.
    If there are multiple edges with the same lowest-weight, it checks that
    none of the vertices are repeated (to maintain matching conditions).
    Since it is a greedy solution, it will at times fail to find the optimal
    lowest-weight matching set.
    '''
    #print "---------------- DEBUGGING GREEDYMATCHING -----------------\n"
    shortestEdge = float('inf')
    mGraph = copy.deepcopy(oddGraph)
    verticesUsed = []
    minWeightGraph = {}

    while len(mGraph) > 0:
        vertexUsed = ''
        currentV = mGraph.keys()[0]
        adjCities = mGraph.pop(currentV)
        nearestNeighbor = 0

        for neighbor in adjCities:
            if neighbor not in verticesUsed and neighbor in vwod:
                if adjCities[neighbor] < shortestEdge:
                    shortestEdge = adjCities[neighbor]
                    nearestNeighbor = neighbor

        verticesUsed.append(nearestNeighbor)
        verticesUsed.append(currentV)
        minWeightGraph[nearestNeighbor] = {currentV: shortestEdge}
        minWeightGraph[currentV] = {nearestNeighbor: shortestEdge}
        mGraph.pop(nearestNeighbor)

        shortestEdge = float('inf')
    #print "MINWEIGHTGRAPH VALUE: {}".format(minWeightGraph)
    # for i in vwod:
        # if i not in minWeightGraph:
        #     print str(i) + " is not in minWeightGraph."
        # for j in minWeightGraph[i]:
        #     if j not in minWeightGraph:
        #         print str(j) + " is adj to " + str(i) + ", but not the other way."

    #print "\n-------------- END DEBUGGING GREEDYMATCHING ---------------\n\n"

    return minWeightGraph

def combine(MSTree, matching):
    multiGraph = copy.deepcopy(MSTree)

    for city in multiGraph:
        for neighbor in multiGraph[city]:
            multiGraph[city][neighbor] = 1
    for city in matching:
        for neighbor in matching[city]:
            if neighbor not in multiGraph[city]:
                multiGraph[city][neighbor] = 1
            else:
                multiGraph[city][neighbor] += 1
    return multiGraph

'''
Takes in a list of vertices representing an Euler Tour.
Returns a list of a tour with no recurring vertices.
Resulting list represents solution to TSP problem.
'''
def makeTSPList(eulerList):
    TSPList = []
    for v in eulerList:
        if v not in TSPList:
            TSPList.append(v)
    return TSPList

def getTSPTourLength(originalGraph, TSPList):
    assert(len(TSPList) == len(originalGraph))
    totalDist = 0
    for i in range(len(TSPList)-1): # i goes from 0 to the element before last in TSPList
        currCity = TSPList[i] #Goes from index 0 to second-to-last
        nextCity = TSPList[i+1] #Goes from index 1 to last
        totalDist += originalGraph[currCity][nextCity]
    firstCity = TSPList[0]
    lastCity = TSPList[-1]
    totalDist += originalGraph[firstCity][lastCity]
    return totalDist


def christofidesTSP(cities, inputFilename):

    #Builds a complete graph with all cities connected
    initialGraph = buildCompleteGraph(cities)

    #Keep original graph, use copy
    copyOfInitialGraph = copy.deepcopy(initialGraph)

    #Get a graph with MST data (predecessor and distanceTo)
    postPrimGraph = prims_mst(copyOfInitialGraph, 0) #MIGHT MAKE THIS RANDOM AND TRY VARIOUS CITIES TO FIND BEST START POINT

    #Extract just the MS tree from the graph
    MSTree = getMSTree(postPrimGraph)

    #Get a list of vertices with odd degree from the MSTree
    vwod = getVerticesWithOddDegree(MSTree)

    #Reduce the initial graph to include only nodes with odd degree from MS Tree
    reducedGraph = reduceG(initialGraph, vwod)

    #Calculate Matching goes here
    matching = greedyMatching(reducedGraph, vwod)

    #Unite Matching and MSTree goes here
    multiGraph = combine(MSTree, matching)

    eTour = []
    #Do Euler Tour on union of Matching and MS Tree
    eulerTour(multiGraph, eTour, 0)

    #Make Euler Circuit Hamiltonian
    TSPList = makeTSPList(eTour)

    #Get TS tour distance
    tourLength = getTSPTourLength(initialGraph, TSPList)

    #Create output file
    createOutputFile(TSPList, tourLength, inputFilename)

    return [TSPList, tourLength]



#--------------------------- SCRIPT STARTS HERE ----------------------------------------------------

inputFilename = sys.argv[1]
cities = getInputData(inputFilename)
totalTime = 0

#For large data sets, use nearest neighbor
if len(cities) > 500: # change to maximum for Christofides
    i = 0
    start = time.clock()
    TourArray = nearestNeighbor(cities, i)
    totalTime = time.clock() - start
    Tour = TourArray[0]
    TourDist = TourArray[1]
    totalTime = time.clock() - start
    createOutputFile(Tour, TourDist, inputFilename)

else:
    start = time.clock()
    TSPTour = christofidesTSP(cities, inputFilename)
    totalTime = time.clock() - start
print "Time: " + str(totalTime) + " sec"

