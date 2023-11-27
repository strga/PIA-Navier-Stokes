#include <fstream>
#include <iostream>
#include "write_to_file.h"

using namespace std;

void writeDataToFile(const vector<vector<double>>& data, const string& filename) {
    ofstream file(filename);
    if (!file.is_open()) {
        cerr << "Failed to open file '" << filename << "' for writing." << endl;
        return;
    }

    for (const auto& row : data) {
        for (size_t j = 0; j < row.size(); ++j) {
            file << row[j];
            if (j < row.size() - 1) file << ",";
        }
        file << "\n";
    }

    file.close();

    if (file.good()) {
        cout << "Data successfully written to '" << filename << "'." << endl;
    } else {
        cerr << "Error occurred when writing to file '" << filename << "'." << endl;
    }
}