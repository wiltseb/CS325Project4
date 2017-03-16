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

# def findExistingMatch():
#     existingMatch = -1
#
#     return existingMatch
#
# def greedyMatching(oddGraph, vwod):
#     '''
#     Greedy algorithm that finds the lowest-weight matching edges.
#     If there are multiple edges with the same lowest-weight, it checks that
#     none of the vertices are repeated (to maintain matching conditions).
#     Since it is a greedy solution, it will at times fail to find the optimal
#     lowest-weight matching set.
#     '''
#     print "---------------- DEBUGGING GREEDYMATCHING -----------------\n"
#     shortest = -1
#     mGraph = oddGraph.copy()
#     verticesUsed = []
#     minWeightGraph = {}
#
#     while len(mGraph) > 0:
#         pairFound = False
#         vertexUsed = ''
#         currentV = mGraph.keys()[0]
#         paths = mGraph.pop(currentV)
#         matchFound = False
#         prevMatch = findExistingMatch(minWeightGraph, currentV)
#         if prevMatch != -1:
#             print "Previous match is {}!".format(prevMatch)
#             minWeightGraph[currentV] = {prevMatch: paths[prevMatch]}
#             vertexUsed = vertex
#             print "Adding {} to {} with weight {} to graph!".format(currentV, prevMatch, paths[prevMatch])
#             verticesUsed.append(vertexUsed)
#             matchFound = True
#             if vertex == 1 and prevMatch == 0:
#                 #print
#                 print currentV, paths[vertex], minWeightGraph[currentV]
#
#         else:
#             for vertex in paths:
#
#                     #continue
#                     #print "No previous match!"
#                     #break
#                 #else:
#                     #x = 1
#
#                 if vertex in vwod and vertex not in verticesUsed:
#                     #print "{} to {}: {}".format(currentV, vertex, paths[vertex])
#                     if vertex in minWeightGraph:
#                         if currentV in minWeightGraph[vertex]:
#                             edgeToAdd = {vertex: paths[vertex]}
#                             minWeightGraph[currentV] = edgeToAdd
#                             shortest = paths[vertex]
#                             vertexUsed = vertex
#
#                     elif shortest == -1:
#                         edgeToAdd = {vertex: paths[vertex]}
#                         minWeightGraph[currentV] = edgeToAdd
#                         shortest = paths[vertex]
#                         vertexUsed = vertex
#
#                     elif paths[vertex] < shortest:
#                         edgeToAdd = {vertex: paths[vertex]}
#                         minWeightGraph[currentV] = edgeToAdd
#                         shortest = paths[vertex]
#                         vertexUsed = vertex
#
#                     verticesUsed.append(vertexUsed)
#             shortest = -1
#     print "MINWEIGHTGRAPH VALUE: {}".format(minWeightGraph)
#     print "\n-------------- END DEBUGGING GREEDYMATCHING ---------------\n\n"
#     return minWeightGraph

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
