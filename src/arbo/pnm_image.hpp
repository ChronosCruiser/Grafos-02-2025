#ifndef PNM_IMAGE_HPP
#define PNM_IMAGE_HPP

#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <cstdint>
#include <iostream>
#include <cmath>

struct Pixel {
    uint8_t r, g, b;
    
    Pixel() : r(0), g(0), b(0) {}
    Pixel(uint8_t r, uint8_t g, uint8_t b) : r(r), g(g), b(b) {}
    
    // Calculate color difference between two pixels
    double colorDifference(const Pixel& other) const {
        double dr = static_cast<double>(r) - static_cast<double>(other.r);
        double dg = static_cast<double>(g) - static_cast<double>(other.g);
        double db = static_cast<double>(b) - static_cast<double>(other.b);
        return std::sqrt(dr * dr + dg * dg + db * db);
    }
};

class PNMImage {
public:
    enum class Format { P3, P6, P2, P5 }; // PPM ASCII, PPM Binary, PGM ASCII, PGM Binary
    
private:
    int width_;
    int height_;
    int maxVal_;
    Format format_;
    std::vector<Pixel> pixels_;
    
    void skipComments(std::ifstream& file) {
        char c;
        while (file.peek() == '#' || std::isspace(file.peek())) {
            if (file.peek() == '#') {
                std::string line;
                std::getline(file, line);
            } else {
                file.get(c);
            }
        }
    }
    
public:
    PNMImage() : width_(0), height_(0), maxVal_(255), format_(Format::P3) {}
    
    bool load(const std::string& filename) {
        std::ifstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot open file " << filename << std::endl;
            return false;
        }
        
        std::string magicNumber;
        file >> magicNumber;
        
        if (magicNumber == "P3") {
            format_ = Format::P3;
        } else if (magicNumber == "P6") {
            format_ = Format::P6;
        } else if (magicNumber == "P2") {
            format_ = Format::P2;
        } else if (magicNumber == "P5") {
            format_ = Format::P5;
        } else {
            std::cerr << "Error: Unsupported PNM format: " << magicNumber << std::endl;
            return false;
        }
        
        skipComments(file);
        file >> width_;
        skipComments(file);
        file >> height_;
        skipComments(file);
        file >> maxVal_;
        
        // Skip single whitespace after maxVal
        file.get();
        
        pixels_.resize(width_ * height_);
        
        if (format_ == Format::P3) {
            // ASCII PPM
            for (int i = 0; i < width_ * height_; ++i) {
                int r, g, b;
                file >> r >> g >> b;
                pixels_[i] = Pixel(
                    static_cast<uint8_t>(r * 255 / maxVal_),
                    static_cast<uint8_t>(g * 255 / maxVal_),
                    static_cast<uint8_t>(b * 255 / maxVal_)
                );
            }
        } else if (format_ == Format::P6) {
            // Binary PPM
            for (int i = 0; i < width_ * height_; ++i) {
                uint8_t rgb[3];
                file.read(reinterpret_cast<char*>(rgb), 3);
                pixels_[i] = Pixel(rgb[0], rgb[1], rgb[2]);
            }
        } else if (format_ == Format::P2) {
            // ASCII PGM
            for (int i = 0; i < width_ * height_; ++i) {
                int gray;
                file >> gray;
                uint8_t g = static_cast<uint8_t>(gray * 255 / maxVal_);
                pixels_[i] = Pixel(g, g, g);
            }
        } else if (format_ == Format::P5) {
            // Binary PGM
            for (int i = 0; i < width_ * height_; ++i) {
                uint8_t gray;
                file.read(reinterpret_cast<char*>(&gray), 1);
                pixels_[i] = Pixel(gray, gray, gray);
            }
        }
        
        std::cout << "Loaded image: " << width_ << "x" << height_ 
                  << " (" << width_ * height_ << " pixels)" << std::endl;
        
        return true;
    }
    
    bool save(const std::string& filename) const {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Error: Cannot create file " << filename << std::endl;
            return false;
        }
        
        // Save as P6 (binary PPM)
        file << "P6\n";
        file << width_ << " " << height_ << "\n";
        file << "255\n";
        
        for (const auto& pixel : pixels_) {
            file.put(pixel.r);
            file.put(pixel.g);
            file.put(pixel.b);
        }
        
        return true;
    }
    
    int getWidth() const { return width_; }
    int getHeight() const { return height_; }
    int getPixelCount() const { return width_ * height_; }
    
    const Pixel& getPixel(int x, int y) const {
        return pixels_[y * width_ + x];
    }
    
    Pixel& getPixel(int x, int y) {
        return pixels_[y * width_ + x];
    }
    
    const Pixel& getPixelByIndex(int index) const {
        return pixels_[index];
    }
    
    Pixel& getPixelByIndex(int index) {
        return pixels_[index];
    }
    
    int pixelIndex(int x, int y) const {
        return y * width_ + x;
    }
    
    std::pair<int, int> indexToCoord(int index) const {
        return {index % width_, index / width_};
    }
};

#endif // PNM_IMAGE_HPP