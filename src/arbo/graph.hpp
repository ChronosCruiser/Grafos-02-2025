#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <vector>
#include <limits>
#include <iostream>
#include "pnm_image.hpp"

struct Edge {
    int from;
    int to;
    double weight;
    
    Edge(int f, int t, double w) : from(f), to(t), weight(w) {}
};

class DirectedWeightedGraph {
private:
    int numVertices_;
    std::vector<Edge> edges_;
    std::vector<std::vector<int>> adjacencyList_; // indices into edges_
    std::vector<std::vector<int>> incomingEdges_; // indices of edges coming into each vertex
    
public:
    DirectedWeightedGraph(int vertices) 
        : numVertices_(vertices), 
          adjacencyList_(vertices),
          incomingEdges_(vertices) {}
    
    void addEdge(int from, int to, double weight) {
        int edgeIndex = edges_.size();
        edges_.emplace_back(from, to, weight);
        adjacencyList_[from].push_back(edgeIndex);
        incomingEdges_[to].push_back(edgeIndex);
    }
    
    int getNumVertices() const { return numVertices_; }
    int getNumEdges() const { return edges_.size(); }
    
    const std::vector<Edge>& getEdges() const { return edges_; }
    const Edge& getEdge(int index) const { return edges_[index]; }
    
    const std::vector<int>& getOutgoingEdges(int vertex) const {
        return adjacencyList_[vertex];
    }
    
    const std::vector<int>& getIncomingEdges(int vertex) const {
        return incomingEdges_[vertex];
    }
    
    // Create graph from PNM image
    static DirectedWeightedGraph fromImage(const PNMImage& image) {
        int width = image.getWidth();
        int height = image.getHeight();
        int numPixels = width * height;
        
        DirectedWeightedGraph graph(numPixels);
        
        // Direction offsets: right, down, left, up
        const int dx[] = {1, 0, -1, 0};
        const int dy[] = {0, 1, 0, -1};
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int currentIndex = image.pixelIndex(x, y);
                const Pixel& currentPixel = image.getPixel(x, y);
                
                // Check all 4 neighbors
                for (int d = 0; d < 4; ++d) {
                    int nx = x + dx[d];
                    int ny = y + dy[d];
                    
                    // Check bounds
                    if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                        int neighborIndex = image.pixelIndex(nx, ny);
                        const Pixel& neighborPixel = image.getPixel(nx, ny);
                        
                        // Weight is the color difference
                        double weight = currentPixel.colorDifference(neighborPixel);
                        
                        // Add directed edge from current to neighbor
                        graph.addEdge(currentIndex, neighborIndex, weight);
                    }
                }
            }
        }
        
        std::cout << "Created graph with " << graph.getNumVertices() 
                  << " vertices and " << graph.getNumEdges() << " edges" << std::endl;
        
        return graph;
    }
};

#endif // GRAPH_HPP