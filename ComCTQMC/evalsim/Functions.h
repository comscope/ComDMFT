#ifndef FUNCTIONS
#define FUNCTIONS

#include "../include/measurements/Measurements.h"
#include "../include/io/Vector.h"
#include "../include/io/Matrix.h"
#include "../include/atomic/GenerateAtomic.h"


io::cmat get_matrix(std::map<std::string, std::complex<double>> const& entries, jsx::value const& jMatrix)
{
    io::cmat matrix(jMatrix.size(), jMatrix.size());
    
    for(std::size_t i = 0; i < jMatrix.size(); ++i)
        for(std::size_t j = 0; j < jMatrix.size(); ++j) {
            auto const entry = jMatrix(i)(j).string();
            
            if(entry != "") matrix(i, j) = entries.at(entry);
        }
    
    //make symmetric if real
    
    return matrix;
}

inline std::vector<io::cmat> get_function_matrix(std::map<std::string, io::cvec> const& functions, jsx::value const& jMatrix)
{
    std::size_t size = std::numeric_limits<std::size_t>::max();
    
    for(auto const& function : functions)
        size = std::min(size, function.second.size());

    std::vector<io::cmat> functionMatrix;
    for(std::size_t n = 0; n < size; ++n) {
        std::map<std::string, std::complex<double>> entries;
        
        for(auto& function : functions)
            entries[function.first] = function.second[n];
        
        functionMatrix.push_back(get_matrix(entries, jMatrix));
    }
    
    return functionMatrix;
}


inline std::map<std::string, std::complex<double>> get_entries(io::cmat const& matrix, jsx::value const& jMatrix)
{
    std::map<std::string, std::pair<int, std::complex<double>>> temp;
    for(std::size_t i = 0; i < jMatrix.size(); ++i)
        for(std::size_t j = 0; j < jMatrix.size(); ++j) {
            auto const entry = jMatrix(i)(j).string();
            
            if(entry != "") {
                if(!temp.count(entry))
                    temp[entry] = std::pair<int, std::complex<double>>(0, .0);
                
                temp.at(entry).first  += 1;
                temp.at(entry).second += matrix(i, j);
            }
        }
    
    std::map<std::string, std::complex<double>> entries;
    for(auto const& entry : temp)
        entries[entry.first] = entry.second.second/std::complex<double>(entry.second.first);
    
    return entries;
}

inline std::map<std::string, io::cvec> get_function_entries(std::vector<io::cmat> const& functionMatrix, jsx::value const& jMatrix)
{
    std::map<std::string, io::cvec> functionEntries;
    
    for(auto const& matrix : functionMatrix) {
        std::map<std::string, std::complex<double>> entries = get_entries(matrix, jMatrix);
        for(auto const& entry : entries)
            functionEntries[entry.first].push_back(entry.second);
    }
    
    return functionEntries;
}

#endif //EVALSIM










