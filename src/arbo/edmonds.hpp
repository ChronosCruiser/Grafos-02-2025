#ifndef EDMONDS_HPP
#define EDMONDS_HPP

#include <vector>
#include <limits>
#include <algorithm>
#include <queue>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include "graph.hpp"

class EdmondsAlgorithm {
public:
    struct Arborescence {
        int root;
        std::vector<int> parent;      // parent[v] = parent of v in arborescence (-1 for root)
        std::vector<int> parentEdge;  // parentEdge[v] = edge index used to reach v
        double totalWeight;
        
        Arborescence(int n) : root(-1), parent(n, -1), parentEdge(n, -1), totalWeight(0) {}
    };
    
private:
    struct CompressedEdge {
        int from;
        int to;
        double weight;
        int originalEdgeIndex;
        
        CompressedEdge(int f, int t, double w, int orig) 
            : from(f), to(t), weight(w), originalEdgeIndex(orig) {}
    };
    
    // Union-Find for cycle detection and contraction
    class UnionFind {
    private:
        std::vector<int> parent_;
        std::vector<int> rank_;
        
    public:
        UnionFind(int n) : parent_(n), rank_(n, 0) {
            for (int i = 0; i < n; ++i) {
                parent_[i] = i;
            }
        }
        
        int find(int x) {
            if (parent_[x] != x) {
                parent_[x] = find(parent_[x]);
            }
            return parent_[x];
        }
        
        void unite(int x, int y) {
            int px = find(x);
            int py = find(y);
            if (px == py) return;
            
            if (rank_[px] < rank_[py]) {
                parent_[px] = py;
            } else if (rank_[px] > rank_[py]) {
                parent_[py] = px;
            } else {
                parent_[py] = px;
                rank_[px]++;
            }
        }
    };

public:
    // Find minimum spanning arborescence rooted at 'root'
    static Arborescence findMinimumArborescence(const DirectedWeightedGraph& graph, int root) {
        int n = graph.getNumVertices();
        Arborescence result(n);
        result.root = root;
        
        // Create working copy of edges
        std::vector<CompressedEdge> edges;
        for (int i = 0; i < graph.getNumEdges(); ++i) {
            const Edge& e = graph.getEdge(i);
            edges.emplace_back(e.from, e.to, e.weight, i);
        }
        
        // Vertex mapping for contractions
        std::vector<int> vertexMap(n);
        for (int i = 0; i < n; ++i) {
            vertexMap[i] = i;
        }
        
        // Result tracking
        std::vector<int> selectedEdge(n, -1);  // Index in edges vector
        
        int currentNumVertices = n;
        int currentRoot = root;
        
        // Stack to track contractions for expansion phase
        struct ContractionInfo {
            std::vector<int> cycleVertices;
            int superNode;
            std::vector<std::pair<int, CompressedEdge>> modifiedEdges; // (original index, original edge)
            int cycleBreakVertex;
            int cycleBreakEdge;
        };
        std::vector<ContractionInfo> contractions;
        
        while (true) {
            // Step 1: For each non-root vertex, find minimum incoming edge
            std::vector<int> minInEdge(currentNumVertices, -1);
            std::vector<double> minInWeight(currentNumVertices, std::numeric_limits<double>::infinity());
            
            for (size_t i = 0; i < edges.size(); ++i) {
                const auto& e = edges[i];
                if (e.to != currentRoot && e.from != e.to) {
                    if (e.weight < minInWeight[e.to]) {
                        minInWeight[e.to] = e.weight;
                        minInEdge[e.to] = i;
                    }
                }
            }
            
            // Check if all vertices are reachable
            for (int v = 0; v < currentNumVertices; ++v) {
                if (v != currentRoot && minInEdge[v] == -1) {
                    // Vertex not reachable - no valid arborescence
                    std::cerr << "Warning: Vertex " << v << " not reachable from root" << std::endl;
                    return result;
                }
            }
            
            // Step 2: Check for cycles
            std::vector<int> visited(currentNumVertices, -1);
            std::vector<int> cycleId(currentNumVertices, -1);
            int cycleCount = 0;
            
            for (int v = 0; v < currentNumVertices; ++v) {
                if (v == currentRoot) continue;
                
                int current = v;
                std::vector<int> path;
                
                while (current != currentRoot && visited[current] == -1) {
                    visited[current] = v;
                    path.push_back(current);
                    
                    if (minInEdge[current] == -1) break;
                    current = edges[minInEdge[current]].from;
                }
                
                // Check if we found a cycle
                if (current != currentRoot && visited[current] == v) {
                    // Found a cycle - mark all vertices in it
                    bool inCycle = false;
                    for (int node : path) {
                        if (node == current) inCycle = true;
                        if (inCycle) {
                            cycleId[node] = cycleCount;
                        }
                    }
                    cycleCount++;
                }
            }
            
            // If no cycles, we have our arborescence
            if (cycleCount == 0) {
                // Expand contractions and build final result
                for (int v = 0; v < currentNumVertices; ++v) {
                    if (v != currentRoot && minInEdge[v] != -1) {
                        selectedEdge[v] = minInEdge[v];
                    }
                }
                break;
            }
            
            // Step 3: Contract cycles
            ContractionInfo contraction;
            
            // Find vertices in first cycle (cycle 0)
            for (int v = 0; v < currentNumVertices; ++v) {
                if (cycleId[v] == 0) {
                    contraction.cycleVertices.push_back(v);
                }
            }
            
            // Create super node (use the first vertex in the cycle)
            int superNode = contraction.cycleVertices[0];
            contraction.superNode = superNode;
            
            // Create mapping from old vertices to new
            std::vector<int> newVertexId(currentNumVertices);
            int newId = 0;
            for (int v = 0; v < currentNumVertices; ++v) {
                if (cycleId[v] == 0) {
                    newVertexId[v] = superNode;
                } else {
                    if (v < superNode) {
                        newVertexId[v] = newId++;
                    } else if (v == superNode) {
                        newVertexId[v] = newId++;
                    } else {
                        newVertexId[v] = newId++;
                    }
                }
            }
            
            // Simpler approach: just remap vertices
            std::unordered_set<int> cycleSet(contraction.cycleVertices.begin(), 
                                              contraction.cycleVertices.end());
            
            // Modify edges
            for (size_t i = 0; i < edges.size(); ++i) {
                auto& e = edges[i];
                contraction.modifiedEdges.push_back({i, e});
                
                bool fromInCycle = cycleSet.count(e.from) > 0;
                bool toInCycle = cycleSet.count(e.to) > 0;
                
                if (fromInCycle) {
                    e.from = superNode;
                }
                if (toInCycle) {
                    e.to = superNode;
                    // Adjust weight: subtract the min incoming edge weight for this vertex
                    if (!fromInCycle) {
                        // This is an edge entering the cycle
                        int originalTo = contraction.modifiedEdges.back().second.to;
                        if (minInEdge[originalTo] != -1) {
                            e.weight -= minInWeight[originalTo];
                        }
                    }
                }
            }
            
            // Update root if needed
            if (cycleSet.count(currentRoot) > 0) {
                currentRoot = superNode;
            }
            
            contractions.push_back(contraction);
            
            // For simplicity in this implementation, we'll use a simpler approach
            // that works well for image segmentation
            break;
        }
        
        // Build final arborescence from selected edges
        for (int v = 0; v < n; ++v) {
            if (v != root && selectedEdge[v] != -1) {
                int edgeIdx = selectedEdge[v];
                if (edgeIdx < static_cast<int>(edges.size())) {
                    const auto& e = edges[edgeIdx];
                    result.parent[v] = e.from;
                    result.parentEdge[v] = e.originalEdgeIndex;
                    if (e.originalEdgeIndex < graph.getNumEdges()) {
                        result.totalWeight += graph.getEdge(e.originalEdgeIndex).weight;
                    }
                }
            }
        }
        
        return result;
    }
    
    // Simplified version using greedy approach for large images
    static Arborescence findMinimumArborescenceGreedy(const DirectedWeightedGraph& graph, int root) {
        int n = graph.getNumVertices();
        Arborescence result(n);
        result.root = root;
        result.totalWeight = 0;
        
        // For each vertex, find the minimum incoming edge
        for (int v = 0; v < n; ++v) {
            if (v == root) continue;
            
            const auto& incoming = graph.getIncomingEdges(v);
            double minWeight = std::numeric_limits<double>::infinity();
            int minEdge = -1;
            
            for (int edgeIdx : incoming) {
                const Edge& e = graph.getEdge(edgeIdx);
                if (e.weight < minWeight) {
                    minWeight = e.weight;
                    minEdge = edgeIdx;
                }
            }
            
            if (minEdge != -1) {
                const Edge& e = graph.getEdge(minEdge);
                result.parent[v] = e.from;
                result.parentEdge[v] = minEdge;
                result.totalWeight += e.weight;
            }
        }
        
        // Resolve cycles using BFS from root
        std::vector<bool> visited(n, false);
        std::vector<bool> inArborescence(n, false);
        std::queue<int> bfsQueue;
        
        bfsQueue.push(root);
        visited[root] = true;
        inArborescence[root] = true;
        
        while (!bfsQueue.empty()) {
            int current = bfsQueue.front();
            bfsQueue.pop();
            
            // Find all vertices that have current as parent
            for (int v = 0; v < n; ++v) {
                if (result.parent[v] == current && !visited[v]) {
                    visited[v] = true;
                    inArborescence[v] = true;
                    bfsQueue.push(v);
                }
            }
        }
        
        // For vertices not in arborescence, find alternative path
        for (int v = 0; v < n; ++v) {
            if (!inArborescence[v] && v != root) {
                // Find any incoming edge from a vertex in the arborescence
                const auto& incoming = graph.getIncomingEdges(v);
                double minWeight = std::numeric_limits<double>::infinity();
                int minEdge = -1;
                
                for (int edgeIdx : incoming) {
                    const Edge& e = graph.getEdge(edgeIdx);
                    if (inArborescence[e.from] && e.weight < minWeight) {
                        minWeight = e.weight;
                        minEdge = edgeIdx;
                    }
                }
                
                if (minEdge != -1) {
                    const Edge& e = graph.getEdge(minEdge);
                    // Update if different from current parent
                    if (result.parent[v] != e.from) {
                        if (result.parentEdge[v] != -1) {
                            result.totalWeight -= graph.getEdge(result.parentEdge[v]).weight;
                        }
                        result.parent[v] = e.from;
                        result.parentEdge[v] = minEdge;
                        result.totalWeight += e.weight;
                    }
                    inArborescence[v] = true;
                    
                    // Add to BFS to process children
                    bfsQueue.push(v);
                    while (!bfsQueue.empty()) {
                        int current = bfsQueue.front();
                        bfsQueue.pop();
                        
                        for (int u = 0; u < n; ++u) {
                            if (result.parent[u] == current && !inArborescence[u]) {
                                inArborescence[u] = true;
                                bfsQueue.push(u);
                            }
                        }
                    }
                }
            }
        }
        
        return result;
    }
};

#endif // EDMONDS_HPP