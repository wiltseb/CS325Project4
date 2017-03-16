'''
function PerfectMatching()
  Input: odds (list of odd vertices), G (adjacency list)
  while !odds.empty do
    v <-- odds.popFront()
    length <-- infinite
    for u (element in) odds do
      if weight(u,v) < length then
        length <-- weight(u,v)
        closest <-- u
      end if
    end for
    G.addEdge(closest,u)
    odds.remove(closest)
  end while
end function
'''
vwod = ['a', 'e', 'd', 'h', 'c', 'x']
Graph = {'a': {'c': 3, 'b': 6, 'e': 9},
        'h': {'b': 5, 'g': 3, 'f': 9},
        'e': {'a': 9, 'c': 9, 'g': 18, 'f': 8},
         'd': {'c': 2, 'b': 2, 'f': 8}}
modifiedGraph = {'a': {'c': 3, 'b': 6, 'e': 9, 'v': 3},
                'h': {'b': 5, 'g': 3, 'f': 9, 'd': 2},
                'e': {'a': 9, 'c': 9, 'g': 18, 'f': 8, 'x': 3},
                'd': {'c': 3, 'b': 2, 'f': 8, 'h': 2},
                'x': {'e': 3, 'z': 3, 'a': 3},
                'c': {'f': 1, 'a': 3, 'b': 2}
                }

# def greedyMatching(oddGraph, vwod):
#     '''
#     Greedy algorithm that finds the lowest-weight matching edges.
#     If there are multiple edges with the same lowest-weight, it checks that
#     none of the vertices are repeated (to maintain matching conditions).
#     Since it is a greedy solution, it will at times fail to find the optimal
#     lowest-weight matching set.
#     '''
#     shortest = -1
#     vUsed = []
#     minWeightGraph = {}
#     for vertex in oddGraph:
#         for v, w in oddGraph[vertex].items():
#             if v in vwod:
#                 if shortest == -1:
#                     edgeToAdd = {v: w}
#                     minWeightGraph[vertex] = edgeToAdd
#                     shortest = w
#                     vUsed.append(vertex)
#                     vUsed.append(v)
#                 elif w == shortest and v not in vUsed and vertex not in vUsed:
#                     edgeToAdd = {v: w}
#                     minWeightGraph[vertex] = edgeToAdd
#                     vUsed.append(vertex)
#                     vUsed.append(v)
#                 elif w < shortest:
#                     minWeightGraph.clear()
#                     del vUsed [:]
#                     edgeToAdd = {v: w}
#                     minWeightGraph[vertex] = edgeToAdd
#                     vUsed.append(vertex)
#                     vUsed.append(v)
#                     shortest = w
#     return minWeightGraph


    # #Find perfect matching with minimum weight
    # matching = greedyMatching(GPrime, vwod)
    # print "Perfect matching with minimum weight from MST vertices with odd degree:"
    # for i in matching:
    #     cityOne = "From city {}".format(i)
    #     cityTwo = "to city {}".format(matching[i].keys()[0])
    #     weight = "with weight {}".format(matching[i].values()[0])
    #     print "{} {} {}".format(cityOne, cityTwo, weight)


def greedyMatching(vwod, Graph):
    shortest = -1
    mGraph = Graph.copy()
    verticesUsed = []
    minWeightGraph = {}
    while len(mGraph) > 0:
        vertexUsed = ''
        currentV = mGraph.keys()[0]
        paths = mGraph.pop(currentV)

        for vertex in paths:
            if vertex in vwod and vertex not in verticesUsed:
                print "{} to {}: {}".format(currentV, vertex, paths[vertex])
                if vertex in minWeightGraph:
                    print minWeightGraph[vertex]
                    if currentV in minWeightGraph[vertex]:
                        edgeToAdd = {vertex: paths[vertex]}
                        print "Adding {} to {}: {} to the minWeightGraph".format(currentV, vertex, paths[vertex])
                        minWeightGraph[currentV] = edgeToAdd
                        shortest = paths[vertex]
                        vertexUsed = vertex
                    else:
                        print "{} already in minWeightGraph and not paired to {}".format(vertex, currentV)
                elif shortest == -1:
                    edgeToAdd = {vertex: paths[vertex]}
                    print "Adding {} to the minWeightGraph".format(edgeToAdd)
                    minWeightGraph[currentV] = edgeToAdd
                    shortest = paths[vertex]
                    vertexUsed = vertex
                elif paths[vertex] < shortest:
                    print "Replacing {} with lower weight path".format(minWeightGraph[currentV])
                    edgeToAdd = {vertex: paths[vertex]}
                    print "Adding {} to {}: {} to the minWeightGraph".format(currentV, vertex, paths[vertex])
                    minWeightGraph[currentV] = edgeToAdd
                    shortest = paths[vertex]
                    vertexUsed = vertex

        verticesUsed.append(vertexUsed)
        print verticesUsed
        shortest = -1
    return minWeightGraph


print vwod
print modifiedGraph

minWeightMatching = greedyMatching(vwod, modifiedGraph)
print minWeightMatching
