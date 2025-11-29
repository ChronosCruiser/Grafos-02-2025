#include <iostream>
#include <string>
#include <chrono>
#include "pnm_image.hpp"
#include "graph.hpp"
#include "edmonds.hpp"
#include "segmentation.hpp"

void printUsage(const char* programName) {
    std::cout << "Usage: " << programName << " [options]\n"
              << "Options:\n"
              << "  -i <file>     Input PNM file (default: test.pnm)\n"
              << "  -o <prefix>   Output file prefix (default: output)\n"
              << "  -t <value>    Threshold for edge cutting (default: 30.0)\n"
              << "  -n <count>    Number of segments (alternative to threshold)\n"
              << "  -m <size>     Minimum segment size (default: 50)\n"
              << "  -h            Show this help message\n";
}

int main(int argc, char* argv[]) {
    // Default parameters
    std::string inputFile = "test.pnm";
    std::string outputPrefix = "output";
    double threshold = 30.0;
    int numSegments = -1;  // -1 means use threshold instead
    int minSegmentSize = 50;
    
    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-i" && i + 1 < argc) {
            inputFile = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            outputPrefix = argv[++i];
        } else if (arg == "-t" && i + 1 < argc) {
            threshold = std::stod(argv[++i]);
        } else if (arg == "-n" && i + 1 < argc) {
            numSegments = std::stoi(argv[++i]);
        } else if (arg == "-m" && i + 1 < argc) {
            minSegmentSize = std::stoi(argv[++i]);
        } else if (arg == "-h") {
            printUsage(argv[0]);
            return 0;
        }
    }
    
    std::cout << "=== Image Segmentation using Edmonds' Algorithm ===" << std::endl;
    std::cout << "Input file: " << inputFile << std::endl;
    
    // Load image
    auto startTime = std::chrono::high_resolution_clock::now();
    
    PNMImage image;
    if (!image.load(inputFile)) {
        std::cerr << "Failed to load image: " << inputFile << std::endl;
        return 1;
    }
    
    auto loadTime = std::chrono::high_resolution_clock::now();
    std::cout << "Image loaded in " 
              << std::chrono::duration_cast<std::chrono::milliseconds>(loadTime - startTime).count()
              << " ms" << std::endl;
    
    // Create graph from image
    std::cout << "\nCreating graph from image..." << std::endl;
    DirectedWeightedGraph graph = DirectedWeightedGraph::fromImage(image);
    
    auto graphTime = std::chrono::high_resolution_clock::now();
    std::cout << "Graph created in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(graphTime - loadTime).count()
              << " ms" << std::endl;
    
    // Perform segmentation
    std::cout << "\nPerforming segmentation..." << std::endl;
    ImageSegmentation::SegmentationResult segmentation;
    
    if (numSegments > 0) {
        std::cout << "Using fixed number of segments: " << numSegments << std::endl;
        segmentation = ImageSegmentation::segmentByEdgeCutting(image, graph, numSegments);
    } else {
        std::cout << "Using threshold: " << threshold << std::endl;
        segmentation = ImageSegmentation::segmentByArborescence(image, graph, threshold, minSegmentSize);
    }
    
    auto segmentTime = std::chrono::high_resolution_clock::now();
    std::cout << "Segmentation completed in "
              << std::chrono::duration_cast<std::chrono::milliseconds>(segmentTime - graphTime).count()
              << " ms" << std::endl;
    std::cout << "Number of segments: " << segmentation.numSegments << std::endl;
    
    // Generate output images
    std::cout << "\nGenerating output images..." << std::endl;
    
    // Colored segments
    std::string segmentedFile = outputPrefix + "_segmented.pnm";
    PNMImage segmentedImage = ImageSegmentation::createSegmentationImage(image, segmentation);
    if (segmentedImage.save(segmentedFile)) {
        std::cout << "Saved: " << segmentedFile << std::endl;
    }
    
    // Boundary overlay
    std::string boundaryFile = outputPrefix + "_boundaries.pnm";
    PNMImage boundaryImage = ImageSegmentation::createBoundaryImage(image, segmentation);
    if (boundaryImage.save(boundaryFile)) {
        std::cout << "Saved: " << boundaryFile << std::endl;
    }
    
    // Average color per segment
    std::string averageFile = outputPrefix + "_average.pnm";
    PNMImage averageImage = ImageSegmentation::createAverageColorImage(image, segmentation);
    if (averageImage.save(averageFile)) {
        std::cout << "Saved: " << averageFile << std::endl;
    }
    
    auto endTime = std::chrono::high_resolution_clock::now();
    std::cout << "\nTotal processing time: "
              << std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count()
              << " ms" << std::endl;
    
    // Print segment statistics
    std::cout << "\n=== Segment Statistics ===" << std::endl;
    std::vector<int> segmentSizes(segmentation.numSegments, 0);
    for (int label : segmentation.labels) {
        segmentSizes[label]++;
    }
    
    int minSize = *std::min_element(segmentSizes.begin(), segmentSizes.end());
    int maxSize = *std::max_element(segmentSizes.begin(), segmentSizes.end());
    double avgSize = static_cast<double>(image.getPixelCount()) / segmentation.numSegments;
    
    std::cout << "Minimum segment size: " << minSize << " pixels" << std::endl;
    std::cout << "Maximum segment size: " << maxSize << " pixels" << std::endl;
    std::cout << "Average segment size: " << avgSize << " pixels" << std::endl;
    
    return 0;
}