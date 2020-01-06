#include <memory>
#include <algorithm>

#include "EvalSim.h"
#include "Functions.h"

#include "../include/measurements/ReadDensityMatrix.h"
#include "../include/measurements/ReadFunction.h"
#include "../include/atomic/GenerateAtomic.h"
#include "../include/linalg/Operators.h"
#include "../include/options/Options.h"


int main(int argc, char** argv)
{
    try {
        if(argc != 2)
            throw std::runtime_error("ReadSim: Wrong number of input parameters !");
        
        
        jsx::value jParams = mpi::read(std::string(argv[1]) + ".json");
        jsx::value jMeasurements = meas::read(mpi::read(std::string(argv[1]) + ".meas.json"));
        
        jParams("hloc") = mpi::read("hloc.json");
        
        jParams["operators"] = ga::construct_annihilation_operators(jParams("hloc"));
        for(auto& jOp : jParams("operators").array()) linalg::make_operator_complex(jOp);
        
        if(jParams.is("dyn")) jParams("dyn") = mpi::read(jParams("dyn").string());
        
        double const beta = jParams("beta").real64(), mu = jParams("mu").real64();
        iOmega const iomega(beta); io::rmat const oneBody = jsx::at<io::rmat>(jParams("hloc")("one body"));
        jsx::value const jHybMatrix = jParams("hybridisation")("matrix");
        
        if( ! jParams.is("complex hybridisation")) jParams["complex hybridisation"] = false;
        
        std::cout << "Reading hybridisation function ... " << std::flush;
        
        std::vector<io::cmat> hyb, hybMoments;
        
        {
            jsx::value jFunctions = mpi::read(jParams("hybridisation")("functions").string());
            
            std::map<std::string, io::cvec> functions;
            for(auto& function : jFunctions.object())
                functions[function.first] = jsx::at<io::cvec>(function.second);
            
            hyb = get_function_matrix(functions, jHybMatrix);
            
            hybMoments = { get_hybridisation_moments(functions, jParams, jHybMatrix) };
        }

        std::cout << "Ok" << std::endl;
        
        
        std::cout << "Reading green function ... " << std::flush;
        
        std::vector<io::cmat> green = read_functions(jMeasurements("green"), jParams, jHybMatrix, hyb.size());

        std::cout << "Ok" << std::endl;
        
        
        std::cout << "Calculating self-energy with dyson ... " << std::flush;
        
        std::vector<io::cmat> selfDyson(green.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
        for(std::size_t n = 0; n < green.size(); ++n) {
            io::cmat green_inv = linalg::inv(green[n]);
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                    selfDyson[n](i, j) = (i == j ? iomega(n) + mu : .0) - oneBody(i, j) - hyb[n](i, j) - green_inv(i, j);
        }
        
        std::cout << "OK" << std::endl;
        
        
        std::vector<io::cmat> self(green.size(), io::cmat(jHybMatrix.size(), jHybMatrix.size()));
        if(jParams.is("green bulla") ? jParams("green bulla").boolean() : true) {
            
            std::cout << "Calculating self-energy with bulla ... " << std::flush;
            
            std::vector<io::cmat> bullaR = read_functions(jMeasurements("bullaR"), jParams, jHybMatrix, hyb.size());
            std::vector<io::cmat> bullaL = read_functions(jMeasurements("bullaL"), jParams, jHybMatrix, hyb.size());

            for(std::size_t n = 0; n < green.size(); ++n) {
                io::cmat green_inv = linalg::inv(green[n]);
                
                linalg::gemm('n', 'n', .5, green_inv, bullaL[n], .0, self[n]);
                linalg::gemm('n', 'n', .5, bullaR[n], green_inv, 1., self[n]);
            }
            
            std::cout << "Ok" << std::endl;
            
        } else
            self = std::move(selfDyson);
        
        
        jsx::value jObservables;
        
        jObservables["sign"] = jMeasurements("sign");
        jObservables["scalar"]["k"] = jsx::at<io::rvec>((jMeasurements("scalar")("k")))[0];
        jObservables["expansion histogram"] = jMeasurements("expansion histogram");
        
        std::cout << "Calculating quantum number observables ... " << std::flush;
        
        opt::complete_qn(jParams);
        
        for(auto jqn : jParams("quantum numbers").object()) {
            auto const& sectorProb = jsx::at<io::rvec>(jMeasurements("sector prob"));
            auto jsqn = ga::construct_sector_qn(jParams("hloc"), jqn.second);
            
            double val = .0, valval = .0;
            for(int s = 0; s < jParams("hloc")("eigen values").size(); ++s) {
                val += sectorProb.at(s)*jsqn(s).real64();
                valval += sectorProb.at(s)*jsqn(s).real64()*jsqn(s).real64();
            }
            
            jObservables("scalar")[jqn.first] = val;
            jObservables("scalar")[jqn.first + jqn.first] = valval;
        }
        
        std::cout << "Ok" << std::endl;
        
        
        jsx::value jDensityMatrix = meas::read_density_matrix(jParams, jMeasurements);
        
        
        std::cout << "Calculating observables ... " << std::flush;
        
        opt::complete_observables(jParams);
        
        for(auto& jOb : jParams("observables").object()) {
            jsx::value jTemp = ga::construct_observable(jParams("hloc"), jOb.second);
            linalg::make_operator_complex(jTemp);
            jObservables("scalar")[jOb.first] = linalg::trace(jDensityMatrix, jTemp).real();
        }
        
        std::cout << "Ok" << std::endl;
        
        
        std::cout << "Calculating green moments ... " << std::flush;
        
        std::vector<io::cmat> greenMoments;
        {
            jsx::value jHamiltonian = get_hamiltonian(jParams);

            jsx::value jC = jsx::array_t(jHybMatrix.size());
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i) {
                linalg::mult('n', 'n',  1., jHamiltonian, jParams("operators")(i), .0, jC(i));
                linalg::mult('n', 'n', -1., jParams("operators")(i), jHamiltonian, 1., jC(i));
            }
            
            greenMoments.resize(3, io::cmat(jHybMatrix.size(), jHybMatrix.size()));

            for(std::size_t i = 0; i < jHybMatrix.size(); ++i) {
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    std::string entry = jHybMatrix(i)(j).string();
                    
                    if(!entry.empty()) {
                        jsx::value temp;
                        
                        linalg::mult('n', 'c', 1., jC(i), jParams("operators")(j), .0, temp);
                        linalg::mult('c', 'n', 1., jParams("operators")(j), jC(i), 1., temp);
                        greenMoments[1](i, j) += linalg::trace(jDensityMatrix, temp);
                        
                        linalg::mult('n', 'c', 1., jC(i), jC(j), .0, temp);
                        linalg::mult('c', 'n', 1., jC(j), jC(i), 1., temp);
                        greenMoments[2](i, j) += linalg::trace(jDensityMatrix, temp);
                    }
                }
                
                greenMoments[0](i, i)  = 1.;
            }

            if(jParams.is("dyn")) {
                double Q1 = -jParams("dyn")(0).real64()*jObservables("scalar")("N").real64(), Q2;

                jsx::value jDensityMatrixPrime = meas::read_density_matrix_prime(jParams, jMeasurements, Q2);

                for(std::size_t i = 0; i < jHybMatrix.size(); ++i) {
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                        std::string entry = jHybMatrix(i)(j).string();
                        
                        if(!entry.empty()) {
                            jsx::value temp;
                            
                            linalg::mult('n', 'c', 1., jC(i), jParams("operators")(j), .0, temp);
                            linalg::mult('c', 'n', 1., jParams("operators")(j), jC(i), 1., temp);
                            
                            linalg::mult('n', 'c', 1., jParams("operators")(i), jC(j), 1., temp);
                            linalg::mult('c', 'n', 1., jC(j), jParams("operators")(i), 1., temp);
                            
                            greenMoments[2](i, j) += linalg::trace(jDensityMatrixPrime, temp);
                        }
                    }
                    
                    greenMoments[1](i, i) += Q1;
                    greenMoments[2](i, i) += Q2;
                }
            }

            greenMoments = get_function_matrix(get_function_entries(greenMoments, jHybMatrix), jHybMatrix);
        }
        
        std::cout << "OK" << std::endl;
        
        
        std::cout << "Calculating self-energy moments ... " << std::flush;
        
        std::vector<io::cmat> selfMoments;
        {
            selfMoments.resize(2, io::cmat(jHybMatrix.size(), jHybMatrix.size()));
            
            io::cmat gm1gm1(jHybMatrix.size(), jHybMatrix.size());
            linalg::gemm('n', 'n', 1., greenMoments[1], greenMoments[1], .0, gm1gm1);
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    selfMoments[0](i, j) += (i == j ? mu : .0) - oneBody(i, j) - greenMoments[1](i, j);
                    selfMoments[1](i, j) += greenMoments[2](i, j) - gm1gm1(i, j);
                }
            
            selfMoments = get_function_matrix(get_function_entries(selfMoments, jHybMatrix), jHybMatrix);
            
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                    greenMoments[2](i, j) += hybMoments[0](i, j);
        }
        
        std::cout << "OK" << std::endl;
        
        
        std::cout << "Adding self-energy high frequency tail ... "  << std::flush;
        
        {
            
            add_self_tail(jHybMatrix, iomega, self, selfMoments, hyb.size());
            
            jObservables["self-energy"] =  write_functions(jParams, jHybMatrix, self, selfMoments);
            
            if(selfDyson.size()) jObservables["self-energy-dyson"] = write_functions(jParams, jHybMatrix, selfDyson, selfMoments);
        }
        
        std::cout << "Ok" << std::endl;
        
        
        std::cout << "Adding green function high frequency tail ... " << std::flush;
        
        {
            for(std::size_t n = green.size(); n < hyb.size(); ++n) {
                io::cmat green_inv(jHybMatrix.size(), jHybMatrix.size());
                
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j)
                        green_inv(i, j) = (i == j ? iomega(n) + mu : .0) - oneBody(i, j) - hyb[n](i, j) - self[n](i, j);

                green.push_back(linalg::inv(green_inv));
            }

            jObservables["green"] = write_functions(jParams, jHybMatrix, green, greenMoments);
        }
        
        std::cout << "Ok" << std::endl;
        
        
        std::cout << "Calculating energy ... " << std::flush;
        
        {
            jsx::value jHamiltonianEff = get_effective_hamiltonian(jParams);
            
            double energy = linalg::trace(jDensityMatrix, jHamiltonianEff).real();
            if(jParams.is("dyn")) energy += jsx::at<io::rvec>(jMeasurements("dynE"))[0];
            jObservables["scalar"]["energy"] = energy;
        }
        
        std::cout << "Ok" << std::endl;
        
        
        std::cout << "Calculating occupation ... " << std::flush;
        
        io::cmat occupation(jHybMatrix.size(), jHybMatrix.size()), correlation(jHybMatrix.size(), jHybMatrix.size());
        {
            jsx::value jn = jsx::array_t(jHybMatrix.size());
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                linalg::mult('c', 'n', 1., jParams("operators")(i), jParams("operators")(i), .0, jn(i));
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    std::string entry = jHybMatrix(i)(j).string();
                    
                    if(!entry.empty()) {
                        jsx::value jOccupation;
                        linalg::mult('c', 'n', 1., jParams("operators")(i), jParams("operators")(j), .0, jOccupation);
                        occupation(i, j) = linalg::trace(jDensityMatrix, jOccupation);
                    }
                    
                    jsx::value jCorrelation;
                    linalg::mult('n', 'n', 1., jn(i), jn(j), .0, jCorrelation);
                    correlation(i, j) = linalg::trace(jDensityMatrix, jCorrelation);
                }
            
            auto occupation_entries = get_entries(occupation, jHybMatrix);
            
            for(auto const& entry : occupation_entries)
                if(jParams("complex hybridisation").boolean()) {
                    jObservables["occupation"][entry.first] = jsx::object_t{{"real", entry.second.real()}, {"imag", entry.second.imag()}};
                } else {
                    jObservables["occupation"][entry.first] = entry.second.real();
                }
            
    
            
            occupation = get_matrix(occupation_entries, jHybMatrix);
        }
        
        std::cout << "Ok" << std::endl;
        
        
        //----------------------------------------------------------------------------------------------------------------------------------------------------
        //----------------------------------------------------------------------------------------------------------------------------------------------------
        //----------------------------------------------------------------------------------------------------------------------------------------------------
        
        
        if(jParams.is("quantum number susceptibility") ? jParams("quantum number susceptibility").boolean() : false) {
            for(auto jqn : jParams("quantum numbers").object()) {
                
                std::cout << "Reading " << jqn.first << " susceptibility ... " << std::flush;
                
                auto const qn = jsx::at<io::rvec>(jqn.second);
                
                double moment = 0;
                io::rvec function(jsx::at<io::rvec>(jMeasurements("susceptibility flavor")("0_0")).size(), .0);
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                        double const fact = qn.at(i)*qn.at(j);
                        auto const meas = jsx::at<io::rvec>(jMeasurements("susceptibility flavor")(std::to_string(i) + "_" + std::to_string(j)));
                        
                        function[0] += fact*meas[0]/(2.*beta);
                        
                        for(std::size_t n = 1; n < function.size(); ++n) function[n] += fact*beta*meas[n]/(4*M_PI*M_PI*n*n);
                        
                        if(i == j) moment -= fact*(jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;
                    }
                function[0] += beta*jObservables("scalar")(jqn.first + jqn.first).real64();
                function[0] -= beta*jObservables("scalar")(jqn.first).real64()*jObservables("scalar")(jqn.first).real64();
                
                std::cout << "Ok" << std::endl;
                
                
                std::cout << "Adding " << jqn.first << " susceptibility high frequency tail ... " << std::flush;
                
                std::size_t const nFit = std::max(static_cast<int>(function.size()/10.), 1);
                
                double A = .0, B = .0;
                for(std::size_t n = function.size() - nFit; n < function.size(); ++n) {
                    A += moment*function[n] + (4*M_PI*M_PI*n*n/(beta*beta))*function[n]*function[n]; B += function[n]*function[n];
                }
                double const alpha = -A/B;
                
                std::size_t const nTail = std::max(static_cast<int>(beta*jParams("susceptibility tail").real64()/(2*M_PI)), 1);
                for(std::size_t n = function.size(); n < nTail; ++n)
                    function.push_back(-moment/((4*M_PI*M_PI*n*n/(beta*beta)) + alpha));
                
                jObservables["susceptibility"][jqn.first]["function"] = function;
                jObservables["susceptibility"][jqn.first]["moment"] = jsx::array_t{{ moment }};
                
                std::cout << "Ok" << std::endl;
            }
        }
    
        
        io::rmat moments(jHybMatrix.size(), jHybMatrix.size());
        if((jParams.is("occupation susceptibility bulla") ? jParams("occupation susceptibility bulla").boolean() : false) ||
           (jParams.is("occupation susceptibility direct") ? jParams("occupation susceptibility direct").boolean() : false))
        {
            std::cout << "Calculating occupation susceptibility moments ... " << std::flush;
            
            jsx::value jHamiltonian = get_hamiltonian(jParams);
            {
                jsx::value jn = jsx::array_t(jHybMatrix.size()), jC = jsx::array_t(jHybMatrix.size());
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i) {
                    linalg::mult('c', 'n', 1., jParams("operators")(i), jParams("operators")(i), .0, jn(i));
                    
                    linalg::mult('n', 'n',  1., jHamiltonian, jn(i), .0, jC(i));
                    linalg::mult('n', 'n', -1., jn(i), jHamiltonian, 1., jC(i));
                }
                
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                        jsx::value moment;
                        
                        linalg::mult('n', 'n',  1., jC(i), jn(j), .0, moment);
                        linalg::mult('n', 'n', -1., jn(j), jC(i), 1., moment);
                        
                        moments(i, j) = linalg::trace(jDensityMatrix, moment).real();
                        if(i == j) moments(i, j) -= (jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;
                    }
            }
            
            std::cout << "Ok" << std::endl;
        }
        
        
        if(jParams.is("occupation susceptibility bulla") ? jParams("occupation susceptibility bulla").boolean() : false) {
            
            std::cout << "Reading bulla occupation susceptibility ... " << std::flush;

            std::vector<std::vector<io::rvec>> susceptibilities(jHybMatrix.size(), std::vector<io::rvec>(jHybMatrix.size()));
            {
                auto constants = moments;
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                        io::rvec measA = jsx::at<io::rvec>(jMeasurements("susceptibility flavor")(std::to_string(i) + "_" + std::to_string(j)));
                        io::rvec measB = jsx::at<io::rvec>(jMeasurements("susceptibility bulla")(std::to_string(i) + "_" + std::to_string(j)));
                        
                        io::rvec function(measA.size(), .0);
                        
                        function[0] = beta*(correlation(i, j).real() - occupation(i, i).real()*occupation(j, j).real()) + (measB[0] + measA[0]/beta)/2.;
                        
                        if(i == j) constants(i, j) += (jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i] + jsx::at<io::rvec>(jMeasurements("flavor k"))[2*i + 1])/beta;

                        for(std::size_t n = 1; n < function.size(); ++n)
                            function[n] = -beta*beta*(constants(i, j) - measB[n] - measA[n]/beta)/(4*M_PI*M_PI*n*n);
                        
                        susceptibilities[i][j] = std::move(function);
                    }
            }

            if( ! jParams("complex hybridisation").boolean()) {
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = i; j < jHybMatrix.size(); ++j)
                        for(std::size_t n = 0; n < susceptibilities[i][j].size(); ++n) {
                            double temp = (susceptibilities[i][j][n] + susceptibilities[j][i][n])/2.;
                            susceptibilities[i][j][n] = susceptibilities[j][i][n] = temp;
                        }
            }
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Adding occupation susceptibility high frequency tail ... " << std::flush;
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j){
                    double const moment = moments(i, j);
                    auto& function = susceptibilities[i][j];
                    std::size_t const nFit = std::max(static_cast<int>(function.size()/10.), 1);
                    
                    double A = .0, B = .0;
                    for(std::size_t n = function.size() - nFit; n < function.size(); ++n) {
                        A += moment*function[n] + (4*M_PI*M_PI*n*n/(beta*beta))*function[n]*function[n]; B += function[n]*function[n];
                    }
                    double const alpha = -A/B;
                    
                    std::size_t const nTail = std::max(static_cast<int>(beta*jParams("susceptibility tail").real64()/(2*M_PI)), 1);
                    for(std::size_t n = function.size(); n < nTail; ++n)
                        function.push_back(-moment/((4*M_PI*M_PI*n*n/(beta*beta)) + alpha));
                }
            
            std::cout << "Ok" << std::endl;
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    jObservables["occupation-susceptibility-bulla"][std::to_string(i) + "_" + std::to_string(j)]["function"] = susceptibilities[i][j];
                    jObservables["occupation-susceptibility-bulla"][std::to_string(i) + "_" + std::to_string(j)]["moment"] = jsx::array_t{{ moments(i, j) }};
                }
        }
        
        
        if(jParams.is("occupation susceptibility direct") ? jParams("occupation susceptibility direct").boolean() : false) {
            
            std::cout << "Reading direct occupation susceptibility ... " << std::flush;
            
            std::vector<std::vector<io::rvec>> susceptibilities(jHybMatrix.size(), std::vector<io::rvec>(jHybMatrix.size()));
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    susceptibilities[i][j] = jsx::at<io::rvec>(jMeasurements("susceptibility direct")(std::to_string(i) + "_" + std::to_string(j)));
                    susceptibilities[i][j][0] -= occupation(i,i).real()*occupation(j, j).real()*beta;
                }
            
            if( ! jParams("complex hybridisation").boolean()) {
                for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                    for(std::size_t j = i; j < jHybMatrix.size(); ++j)
                        for(std::size_t n = 0; n < susceptibilities[i][j].size(); ++n) {
                            double temp = (susceptibilities[i][j][n] + susceptibilities[j][i][n])/2.;
                            susceptibilities[i][j][n] = susceptibilities[j][i][n] = temp;
                        }
            }
            
            for(std::size_t i = 0; i < jHybMatrix.size(); ++i)
                for(std::size_t j = 0; j < jHybMatrix.size(); ++j) {
                    jObservables["occupation-susceptibility-direct"][std::to_string(i) + "_" + std::to_string(j)]["function"] = susceptibilities[i][j];
                    jObservables["occupation-susceptibility-direct"][std::to_string(i) + "_" + std::to_string(j)]["moment"] = jsx::array_t{{ moments(i, j) }};
                }
            
            
            std::cout << "Ok" << std::endl;
        }
        
        
        if(jParams.is("probabilities")) {
            std::cout << "Reading impurity probabilities ... " << std::flush;
            
            jsx::value jHamiltonianEff = get_effective_hamiltonian(jParams);
            
            std::vector<io::rvec> data;
            
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i) {
                    std::vector<double> temp(jParams("probabilities").size() + 1);
                    temp.back() = std::abs(jsx::at<io::cmat>(jDensityMatrix(sector)("matrix"))(i, i).real());     // take abs(real) value ?
                    data.push_back(temp);
                }
            
            for(int qn = 0; qn < jParams("probabilities").size(); ++qn) {
                std::string const name = jParams("probabilities")(qn).string();
                
                if(jParams("quantum numbers").is(name)) {
                    jsx::value jQuantumNumber = ga::construct_sector_qn(jParams("hloc"), jParams("quantum numbers")(name));
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(jQuantumNumber(sector).real64(), 8);
                    
                } else if(jParams("observables").is(name)) {
                    jsx::value jObservable = ga::construct_observable(jParams("hloc"), jParams("observables")(name));
                    
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn]= truncate(jsx::at<io::rmat>(jObservable(sector)("matrix"))(i, i), 8);
                    
                } else if(name == "energy") {
                    int index = 0;
                    for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                        for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i)
                            data[index++][qn] = truncate(jsx::at<io::cmat>(jHamiltonianEff(sector)("matrix"))(i, i).real(), 8);
                    
                } else
                    throw std::runtime_error("quantum number " + name + " for labeling probabilities not found");
            }
            
            if(jParams("probabilities").size()) std::sort(data.begin(), data.end(), VecLess(jParams("probabilities").size() + 1, 0, jParams("probabilities").size()));
            
            jsx::array_t temp;
            for(auto& entry : data) temp.push_back(entry);
            jObservables["probabilities"]["quantum numbers"] = std::move(temp);
            
 
            {
                jsx::value jTransformation = jParams("hloc")("transformation"), jTemp;
                linalg::make_operator_complex(jTransformation);
                
                linalg::mult('n', 'c', 1., jDensityMatrix, jTransformation, .0, jTemp);
                linalg::mult('n', 'n', 1., jTransformation, jTemp, .0, jDensityMatrix);
            }
 
            jsx::value jOccupationStates = ga::construct_occupation_states(jParams("hloc"));
            
            data.clear(); int index = 0;
            for(int sector = 0; sector < jParams("hloc")("eigen values").size(); ++sector)
                for(int i = 0; i < jsx::at<io::rvec>(jParams("hloc")("eigen values")(sector)).size(); ++i) {
                    io::rvec temp = jsx::at<io::rvec>(jOccupationStates(index++));
                    temp.push_back(jsx::at<io::cmat>(jDensityMatrix(sector)("matrix"))(i, i).real());       // take abs(real) value ?
                    data.push_back(temp);
                }

            std::sort(data.begin(), data.end(), VecLess(jHybMatrix.size() + 1, jHybMatrix.size(), jHybMatrix.size() + 1));
            
            temp.clear();
            for(auto& entry : data) temp.push_back(entry);
            jObservables["probabilities"]["occupation numbers"] = std::move(temp);
            
            std::cout << "Ok" << std::endl;
        }
        
        
        mpi::write(jObservables, std::string(argv[1]) + ".obs.json");
    }
    catch(std::exception& exc) {
        std::cerr << exc.what() << "( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
        return -1;
    }
    catch(...) {
        std::cerr << "Fatal Error: Unknown Exception! ( Thrown from worker " << mpi::rank() << " )" << std::endl;
        
        return -2;
    }
    
    return 0;
}












