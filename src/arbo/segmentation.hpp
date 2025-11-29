#ifndef SEGMENTATION_HPP
#define SEGMENTATION_HPP

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <algorithm>
#include <random>
#include <cmath>
#include "pnm_image.hpp"
#include "graph.hpp"
#include "edmonds.hpp"

class ImageSegmentation {
public:
    struct SegmentationResult {
        std::vector<int> labels;           // Segment label for each pixel
        int numSegments;
        std::vector<Pixel> segmentColors;  // Color for each segment (for visualization)
    };
    
private:
    // Generate distinct colors for visualization
    static std::vector<Pixel> generateDistinctColors(int count) {
        std::vector<Pixel> colors;
        std::mt19937 rng(42); // Fixed seed for reproducibility
        
        for (int i = 0; i < count; ++i) {
            // Use HSV to RGB conversion for distinct colors
            float h = static_cast<float>(i) / count;
            float s = 0.7f + (rng() % 30) / 100.0f;
            float v = 0.8f + (rng() % 20) / 100.0f;
            
            // HSV to RGB
            float c = v * s;
            float x = c * (1 - std::abs(std::fmod(h * 6, 2) - 1));
            float m = v - c;
            
            float r, g, b;
            if (h < 1.0f/6) { r = c; g = x; b = 0; }
            else if (h < 2.0f/6) { r = x; g = c; b = 0; }
            else if (h < 3.0f/6) { r = 0; g = c; b = x; }
            else if (h < 4.0f/6) { r = 0; g = x; b = c; }
            else if (h < 5.0f/6) { r = x; g = 0; b = c; }
            else { r = c; g = 0; b = x; }
            
            colors.emplace_back(
                static_cast<uint8_t>((r + m) * 255),
                static_cast<uint8_t>((g + m) * 255),
                static_cast<uint8_t>((b + m) * 255)
            );
        }
        
        return colors;
    }
    
public:
    // Segment image using minimum spanning arborescence
    static SegmentationResult segmentByArborescence(
        const PNMImage& image,
        const DirectedWeightedGraph& graph,
        double threshold,
        int minSegmentSize = 10
    ) {
        int n = graph.getNumVertices();
        SegmentationResult result;
        result.labels.resize(n, -1);
        
        // Find root (use center pixel)
        int centerX = image.getWidth() / 2;
        int centerY = image.getHeight() / 2;
        int root = image.pixelIndex(centerX, centerY);
        
        std::cout << "Finding minimum spanning arborescence from root " << root << "..." << std::endl;
        
        // Find minimum spanning arborescence
        auto arborescence = EdmondsAlgorithm::findMinimumArborescenceGreedy(graph, root);
        
        std::cout << "Arborescence total weight: " << arborescence.totalWeight << std::endl;
        
        // Cut edges with weight above threshold to create segments
        std::vector<int> parent = arborescence.parent;
        
        // Remove edges above threshold
        int edgesCut = 0;
        for (int v = 0; v < n; ++v) {
            if (arborescence.parentEdge[v] != -1) {
                const Edge& e = graph.getEdge(arborescence.parentEdge[v]);
                if (e.weight > threshold) {
                    parent[v] = -1;  // Cut this edge
                    edgesCut++;
                }
            }
        }
        
        std::cout << "Cut " << edgesCut << " edges above threshold " << threshold << std::endl;
        
        // Find connected components (segments)
        int currentLabel = 0;
        for (int v = 0; v < n; ++v) {
            if (result.labels[v] == -1) {
                // BFS to label this component
                std::queue<int> q;
                q.push(v);
                result.labels[v] = currentLabel;
                
                while (!q.empty()) {
                    int current = q.front();
                    q.pop();
                    
                    // Check all neighbors in the arborescence
                    // Children: vertices where parent == current
                    for (int u = 0; u < n; ++u) {
                        if (parent[u] == current && result.labels[u] == -1) {
                            result.labels[u] = currentLabel;
                            q.push(u);
                        }
                    }
                    
                    // Parent
                    if (parent[current] != -1 && result.labels[parent[current]] == -1) {
                        result.labels[parent[current]] = currentLabel;
                        q.push(parent[current]);
                    }
                }
                
                currentLabel++;
            }
        }
        
        result.numSegments = currentLabel;
        std::cout << "Found " << result.numSegments << " initial segments" << std::endl;
        
        // Merge small segments
        if (minSegmentSize > 1) {
            result = mergeSmallSegments(image, graph, result, minSegmentSize);
            std::cout << "After merging: " << result.numSegments << " segments" << std::endl;
        }
        
        // Generate colors for visualization
        result.segmentColors = generateDistinctColors(result.numSegments);
        
        return result;
    }
    
    // Alternative: segment by cutting high-weight edges in arborescence
    static SegmentationResult segmentByEdgeCutting(
        const PNMImage& image,
        const DirectedWeightedGraph& graph,
        int numSegments
    ) {
        int n = graph.getNumVertices();
        SegmentationResult result;
        result.labels.resize(n, -1);
        
        // Find root
        int root = image.pixelIndex(image.getWidth() / 2, image.getHeight() / 2);
        
        // Find minimum spanning arborescence
        auto arborescence = EdmondsAlgorithm::findMinimumArborescenceGreedy(graph, root);
        
        // Collect all edges with their weights
        std::vector<std::pair<double, int>> edgeWeights; // (weight, vertex)
        for (int v = 0; v < n; ++v) {
            if (arborescence.parentEdge[v] != -1) {
                const Edge& e = graph.getEdge(arborescence.parentEdge[v]);
                edgeWeights.push_back({e.weight, v});
            }
        }
        
        // Sort by weight (descending)
        std::sort(edgeWeights.begin(), edgeWeights.end(), std::greater<>());
        
        // Cut the (numSegments - 1) highest weight edges
        std::vector<int> parent = arborescence.parent;
        int toCut = std::min(numSegments - 1, static_cast<int>(edgeWeights.size()));
        
        for (int i = 0; i < toCut; ++i) {
            parent[edgeWeights[i].second] = -1;
        }
        
        // Find connected components
        int currentLabel = 0;
        for (int v = 0; v < n; ++v) {
            if (result.labels[v] == -1) {
                std::queue<int> q;
                q.push(v);
                result.labels[v] = currentLabel;
                
                while (!q.empty()) {
                    int current = q.front();
                    q.pop();
                    
                    for (int u = 0; u < n; ++u) {
                        if (parent[u] == current && result.labels[u] == -1) {
                            result.labels[u] = currentLabel;
                            q.push(u);
                        }
                    }
                    
                    if (parent[current] != -1 && result.labels[parent[current]] == -1) {
                        result.labels[parent[current]] = currentLabel;
                        q.push(parent[current]);
                    }
                }
                
                currentLabel++;
            }
        }
        
        result.numSegments = currentLabel;
        result.segmentColors = generateDistinctColors(result.numSegments);
        
        return result;
    }
    
    // Merge small segments into neighboring segments
    static SegmentationResult mergeSmallSegments(
        const PNMImage& image,
        const DirectedWeightedGraph& graph,
        const SegmentationResult& input,
        int minSize
    ) {
        int n = graph.getNumVertices();
        SegmentationResult result = input;
        
        // Count segment sizes
        std::vector<int> segmentSize(input.numSegments, 0);
        for (int v = 0; v < n; ++v) {
            segmentSize[result.labels[v]]++;
        }
        
        // Create mapping for merged segments
        std::vector<int> segmentMapping(input.numSegments);
        for (int i = 0; i < input.numSegments; ++i) {
            segmentMapping[i] = i;
        }
        
        // Find and merge small segments
        bool changed = true;
        while (changed) {
            changed = false;
            
            // Recalculate sizes
            std::fill(segmentSize.begin(), segmentSize.end(), 0);
            for (int v = 0; v < n; ++v) {
                segmentSize[result.labels[v]]++;
            }
            
            for (int seg = 0; seg < input.numSegments; ++seg) {
                if (segmentSize[seg] > 0 && segmentSize[seg] < minSize) {
                    // Find neighboring segment to merge with
                    std::unordered_map<int, double> neighborWeights;
                    
                    for (int v = 0; v < n; ++v) {
                        if (result.labels[v] != seg) continue;
                        
                        // Check neighbors
                        for (int edgeIdx : graph.getOutgoingEdges(v)) {
                            const Edge& e = graph.getEdge(edgeIdx);
                            int neighborSeg = result.labels[e.to];
                            if (neighborSeg != seg) {
                                neighborWeights[neighborSeg] += 1.0 / (e.weight + 1);
                            }
                        }
                    }
                    
                    // Merge with the most connected neighbor
                    if (!neighborWeights.empty()) {
                        int bestNeighbor = -1;
                        double bestWeight = -1;
                        for (const auto& [neighbor, weight] : neighborWeights) {
                            if (weight > bestWeight && segmentSize[neighbor] > 0) {
                                bestWeight = weight;
                                bestNeighbor = neighbor;
                            }
                        }
                        
                        if (bestNeighbor != -1) {
                            // Merge seg into bestNeighbor
                            for (int v = 0; v < n; ++v) {
                                if (result.labels[v] == seg) {
                                    result.labels[v] = bestNeighbor;
                                }
                            }
                            segmentSize[bestNeighbor] += segmentSize[seg];
                            segmentSize[seg] = 0;
                            changed = true;
                        }
                    }
                }
            }
        }
        
        // Renumber segments to be contiguous
        std::unordered_map<int, int> renumber;
        int newLabel = 0;
        for (int v = 0; v < n; ++v) {
            int oldLabel = result.labels[v];
            if (renumber.find(oldLabel) == renumber.end()) {
                renumber[oldLabel] = newLabel++;
            }
            result.labels[v] = renumber[oldLabel];
        }
        
        result.numSegments = newLabel;
        result.segmentColors = generateDistinctColors(result.numSegments);
        
        return result;
    }
    
    // Create visualization image
    static PNMImage createSegmentationImage(
        const PNMImage& original,
        const SegmentationResult& segmentation
    ) {
        PNMImage result;
        int width = original.getWidth();
        int height = original.getHeight();
        
        // Create a copy
        result = original;
        
        // Color each pixel according to its segment
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = original.pixelIndex(x, y);
                int label = segmentation.labels[idx];
                result.getPixel(x, y) = segmentation.segmentColors[label];
            }
        }
        
        return result;
    }
    
    // Create boundary overlay image
    static PNMImage createBoundaryImage(
        const PNMImage& original,
        const SegmentationResult& segmentation
    ) {
        PNMImage result = original;
        int width = original.getWidth();
        int height = original.getHeight();
        
        // Mark boundary pixels
        Pixel boundaryColor(255, 0, 0); // Red boundaries
        
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = original.pixelIndex(x, y);
                int label = segmentation.labels[idx];
                
                // Check if this is a boundary pixel
                bool isBoundary = false;
                
                const int dx[] = {1, 0, -1, 0};
                const int dy[] = {0, 1, 0, -1};
                
                for (int d = 0; d < 4; ++d) {
                    int nx = x + dx[d];
                    int ny = y + dy[d];
                    
                    if (nx >= 0 && nx < width && ny >= 0 && ny < height) {
                        int neighborIdx = original.pixelIndex(nx, ny);
                        if (segmentation.labels[neighborIdx] != label) {
                            isBoundary = true;
                            break;
                        }
                    }
                }
                
                if (isBoundary) {
                    result.getPixel(x, y) = boundaryColor;
                }
            }
        }
        
        return result;
    }
    
    // Calculate average color for each segment and create smooth visualization
    static PNMImage createAverageColorImage(
        const PNMImage& original,
        const SegmentationResult& segmentation
    ) {
        int width = original.getWidth();
        int height = original.getHeight();
        int n = width * height;
        
        // Calculate average color for each segment
        std::vector<double> sumR(segmentation.numSegments, 0);
        std::vector<double> sumG(segmentation.numSegments, 0);
        std::vector<double> sumB(segmentation.numSegments, 0);
        std::vector<int> count(segmentation.numSegments, 0);
        
        for (int i = 0; i < n; ++i) {
            int label = segmentation.labels[i];
            const Pixel& p = original.getPixelByIndex(i);
            sumR[label] += p.r;
            sumG[label] += p.g;
            sumB[label] += p.b;
            count[label]++;
        }
        
        std::vector<Pixel> avgColors(segmentation.numSegments);
        for (int i = 0; i < segmentation.numSegments; ++i) {
            if (count[i] > 0) {
                avgColors[i] = Pixel(
                    static_cast<uint8_t>(sumR[i] / count[i]),
                    static_cast<uint8_t>(sumG[i] / count[i]),
                    static_cast<uint8_t>(sumB[i] / count[i])
                );
            }
        }
        
        // Create result image
        PNMImage result = original;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int idx = original.pixelIndex(x, y);
                result.getPixel(x, y) = avgColors[segmentation.labels[idx]];
            }
        }
        
        return result;
    }
};

#endif // SEGMENTATION_HPP