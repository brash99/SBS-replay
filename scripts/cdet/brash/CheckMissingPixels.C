#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <sstream>

// Function to get missing pixels for a given layer, submodule, side, and PMT
std::vector<int> GetMissingPixels(int layerNum, int submoduleNum, int sideNum, int pmtNum) {
    std::vector<int> missingPixels;
    std::ifstream file("unusedPixels_parsed.csv");

    if (!file.is_open()) {
        std::cerr << "Error: Could not open unusedPixels_parsed.csv" << std::endl;
        return missingPixels;
    }

    // Convert layer and submodule to module
    int moduleNum = (layerNum == 1) ? submoduleNum : (submoduleNum + 3);

    // Map input side to actual side in file
    int fileSide = (sideNum == 0) ? 2 : 1;

    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        int module, side, pmt, pixel;
        char comma;

        if (ss >> module >> comma >> side >> comma >> pmt >> comma >> pixel) {
            // Match with adjusted side and reversed PMT
            if (module == moduleNum && side == fileSide && pmt == (15 - pmtNum)) {
                if (sideNum == 0) {
                    missingPixels.push_back(pixel);
                } else {
                    missingPixels.push_back(17 - pixel);
                }
            }
        }
    }

    file.close();
    return missingPixels;
}


// Main program
void CheckMissingPixels(bool showAll=true) {
    for (int layer = 1; layer <= 2; ++layer) {
        for (int side = 0; side <= 1; ++side) {
            for (int submodule = 1; submodule <= 3; ++submodule) {
                for (int pmt = 1; pmt <= 14; ++pmt) {
                    // Get the missing pixels for this configuration
                    std::vector<int> missing = GetMissingPixels(layer, submodule, side, pmt);

		    if (showAll){
		      for (int pixel = 1; pixel <= 16; ++pixel) {
                        // Check if this pixel is in the missing list
                        bool isMissing = std::find(missing.begin(), missing.end(), pixel) != missing.end();

                        // Calculate paddle number
                        int paddle = (layer - 1) * 1344 +
                                     (submodule - 1) * 224 +
                                     side * 672 +
                                     (pmt - 1) * 16 +
                                     pixel;

                        // Output result
                        std::cout << "Layer " << layer <<  " Module " << submodule << " Side " << side + 1 << " PMT " << pmt << " pixel " << pixel << " Paddle " << paddle << ": "
                                  << (isMissing ? "MISSING" : "OK") << std::endl;
		      }
		    }
		    else {
		      for (int pixel : missing) {
                        // Calculate paddle number
                        int paddle = (layer - 1) * 1344 +
                                     (submodule - 1) * 224 +
                                     side * 672 +
                                     (pmt - 1) * 16 +
                                     pixel;

                        // Output only missing pixels
                        //std::cout << "Layer " << layer << " Module " << submodule
                        //          << " Side " << side + 1 << " PMT " << pmt
                        //          << " Missing Pixel " << pixel
                        //          << " Paddle " << paddle << std::endl;
			std::cout << paddle - 1 << ", ";
		      }
		    }
                }
            }
        }
    }
}

