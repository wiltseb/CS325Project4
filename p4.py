import math
import sys
import operator
import copy
from collections import defaultdict

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
tour: list city IDs
outFilename: name of output file
'''
def createOutputFile(tourLength, tour, inFilename):
    outFile = inFilename + ".tour"
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
def eulerTour(G):

	eTour = []
	E = copy.deepcopy(G) #added by Brian to avoid changes to one affecting the other

	edgeCount = defaultdict(int)

	def getTour(u):
		for e in E:
			if u == e[0]:
				u,v = e
				E.remove(e)
				getTour(v)

			elif u == e[1]:
				v,u = e
				E.remove(e)
				getTour(v)

		eTour.insert(0,u)

	for x,y in G:               #This threw an error saying 'int' is not iterable; maybe you want a nested loop?
		edgeCount[x] += 1
		edgeCount[y] += 1

	begin = graph[0][0]  #Is this supposed to be G[0][0]? If it is, there's no edge between city 0 and itself.

	for x,y in edgeCount.iteritems():
		if y % 2 > 0:
			begin = x
			break

	curr = begin
	getTour(curr)

	return eTour

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


def findExistingMatch(matches, toMatch):
    '''
    Used to check if a matching pair of vertices already exists
    '''
    existingMatch = -1
    for vertex in matches:
        for v in matches[vertex]:
            if v == toMatch:
                #print "{}".format(vertex)
                existingMatch = vertex
    return existingMatch

def greedyMatching(oddGraph, vwod):
    '''
    Greedy algorithm that finds the lowest-weight matching edges.
    If there are multiple edges with the same lowest-weight, it checks that
    none of the vertices are repeated (to maintain matching conditions).
    Since it is a greedy solution, it will at times fail to find the optimal
    lowest-weight matching set.
    '''
    #print "---------------- DEBUGGING GREEDYMATCHING -----------------\n"
    shortest = -1
    mGraph = oddGraph.copy()
    verticesUsed = []
    minWeightGraph = {}

    while len(mGraph) > 0:
        pairFound = False
        vertexUsed = ''
        currentV = mGraph.keys()[0]
        paths = mGraph.pop(currentV)
        matchFound = False
        prevMatch = findExistingMatch(minWeightGraph, currentV)
        if prevMatch != -1:
            minWeightGraph[currentV] = {prevMatch: paths[prevMatch]}
            if currentV not in verticesUsed:
                verticesUsed.append(currentV)
            if prevMatch not in verticesUsed:
                verticesUsed.append(prevMatch)

        else:
            for vertex in paths:
                if vertex in vwod and vertex not in verticesUsed:
                    if vertex in minWeightGraph:
                        if currentV in minWeightGraph[vertex]:
                            edgeToAdd = {vertex: paths[vertex]}
                            minWeightGraph[currentV] = edgeToAdd
                            shortest = paths[vertex]
                            vertexUsed = vertex

                    elif shortest == -1:
                        edgeToAdd = {vertex: paths[vertex]}
                        minWeightGraph[currentV] = edgeToAdd
                        shortest = paths[vertex]
                        vertexUsed = vertex

                    elif paths[vertex] < shortest:
                        edgeToAdd = {vertex: paths[vertex]}
                        minWeightGraph[currentV] = edgeToAdd
                        shortest = paths[vertex]
                        vertexUsed = vertex

                    verticesUsed.append(vertexUsed)
            shortest = -1
    #print "MINWEIGHTGRAPH VALUE: {}".format(minWeightGraph)
    #print "\n-------------- END DEBUGGING GREEDYMATCHING ---------------\n\n"
    return minWeightGraph

def combine(MSTree, matching):
    multiGraph = copy.deepcopy(MSTree)

    for city in multiGraph:
        for neighbor in multiGraph[city]:
            multiGraph[city][neighbor] = 1
    print multiGraph
    for city in matching:
        for neighbor in matching[city]:
            if neighbor not in multiGraph[city]:
                multiGraph[city][neighbor] = 1
            else:
                multiGraph[city][neighbor] += 1

    print multiGraph
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
    TSPList.append(TSPList[0]) #end at start vertex
    return TSPList

def getTSPTourLength(originalGraph, TSPList):
    assert(len(TSPList) == len(originalGraph))
    totalDist = 0
    for i in range(len(TSPList) - 1): # i goes from 0 to the element before last in TSPList
        currCity = TSPList[i] #Goes from index 0 to second-to-last
        nextCity = TSPList[i+1] #Goes from index 1 to last
        totalDist += originalGraph[currCity][nextCity] #last 'nextCity' is origin
    return totalDist


def solveTSP(inputFilename):
    #Takes in file as described in project specs
    vertices = getInputData(inputFilename)

    #Builds a complete graph with all cities connected
    initialGraph = buildCompleteGraph(vertices)

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

    #Do Euler Tour on union of Matching and MS Tree
    eulerList = eulerTour(multiGraph)

    #Make Euler Circuit Hamiltonian
    TSPList = makeTSPList(eulerList)

    #Get TS tour distance
    tourLength = getTSPTourLength(initialGraph, TSPList)

    #Create output file
    createOutputFile(tourLength, TSPList, outputFilename)




'''
FUNCTION USAGE IN END PRODUCT(Must accept problem instances on command line):
getInputData(sys.arv[1])
TSP function call
createOutputFile(tourLength, tour, sys.argv[2])
'''

'''
#Use to build graph from file data
vertices = getInputData("tsp_example_2.txt");
G = buildCompleteGraph(vertices)

#Input graph to test and city to start MST from
testMSTReduce(testGraph, 'b')
'''
solveTSP("tsp_example_1.txt")
