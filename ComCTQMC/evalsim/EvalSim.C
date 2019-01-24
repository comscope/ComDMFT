#include <memory>

#include "EvalSim.h"


int main(int argc, char** argv)
{
    try {
        if(argc != 2)
            throw std::runtime_error("ReadSim: Wrong number of input parameters !");
        
        std::string name = argv[1];
        jsx::value jSimulation; mpi::read(name + ".meas.json", jSimulation);
        
        jsx::value const& jParams = jSimulation("Parameters");
        jsx::value const jMeasurements = meas::read(jSimulation("Measurements"));
        
        jsx::value jObservables;

        jObservables["sign"] = jMeasurements("Sign");
        jObservables["scalar"] = jMeasurements("Scal");
        jObservables["expansion histogram"] = jMeasurements("ExpansionHist");

        for(auto const& jEntry : jMeasurements("Susc").object()) {
            jObservables["susceptibility"][jEntry.first]["function"] = jEntry.second;
            
            double const val0 = jsx::at<io::rvec>(jObservables("scalar")(jEntry.first))[0];
            jsx::at<io::rvec>(jObservables["susceptibility"][jEntry.first]["function"])[0] -= jParams("beta").real64()*val0*val0;
            
            jObservables["susceptibility"][jEntry.first]["moments"] = std::string("To Do");
        }

        std::vector<io::cmat> green;
        {
            std::map<std::string, io::cvec> functions;

            for(auto& function : jMeasurements("Green").object())
                functions[function.first] = jsx::at<io::cvec>(function.second);

            green = get_function_matrix(functions, jParams("hybridisation")("matrix"));
        }

        iOmega iomega(jParams("beta").real64());
        
        std::vector<io::cmat> hyb;
        std::vector<std::map<std::string, std::complex<double>>> hm_entries(1);
        {
            std::map<std::string, io::cvec> functions;
            
            jsx::value jFunctions; mpi::read(jParams("hybridisation")("functions").string(), jFunctions);
            
            for(auto& function : jFunctions.object())
                functions[function.first] = jsx::at<io::cvec>(function.second);
            
            hyb = get_function_matrix(functions, jParams("hybridisation")("matrix"));
            
            for(auto const& function : functions) { // cross-check with Hyb.h Simple
                std::size_t nFit = std::max(static_cast<int>(function.second.size()*.1), 1);
                
                double ReD = .0, D2 = .0, ReiwD = 0.;
                for(std::size_t n = function.second.size() - nFit; n < function.second.size(); ++n) {
                    ReD += function.second[n].real();
                    D2 += std::abs(function.second[n])*std::abs(function.second[n]);
                    ReiwD += -iomega(n).imag()*function.second[n].imag();
                }
                
                hm_entries[0][function.first] = ReiwD/(nFit - ReD*ReD/D2);
            }
        }

        
        jsx::value jAtomic; std::unique_ptr<Ga::GenerateAtomic> generateAtomic;
        
        if(jParams.is("impurity"))
            mpi::read(jParams("impurity").string(), jAtomic);
        else {
            mpi::read("Hloc.json", jAtomic["hloc"]);
            generateAtomic = std::unique_ptr<Ga::GenerateAtomic>(new Ga::GenerateAtomic(jAtomic["hloc"]));
        }
        

        std::cout << "Calculating self-energy ... " << std::flush;
        
        std::size_t const flavor_number = jParams("hybridisation")("matrix").size();
        
        std::vector<io::cmat> self(std::min(green.size(), hyb.size()));
        for(std::size_t n = 0; n < self.size(); ++n) {
            io::cmat green_inv = linalg::inv(green[n]);
            
            self[n].resize(flavor_number, flavor_number);
            for(std::size_t i = 0; i < flavor_number; ++i)
                for(std::size_t j = 0; j < flavor_number; ++j)
                    self[n](i, j) = (i == j ? iomega(n) + jParams("mu").real64() : .0) - jAtomic("hloc")("one body")(i)(j).real64() - hyb[n](i, j) - green_inv(i, j);
        }
        
        std::cout << "OK" << std::endl;
        
        
        std::vector<std::map<std::string, std::complex<double>>> gm_entries;
        std::vector<std::map<std::string, std::complex<double>>> sm_entries;
        
        if(jMeasurements.is("DensityMatrix")) {
            
            if(generateAtomic.get() != nullptr)
                generateAtomic->write_operators(jAtomic["operators"]);
            
            setup_impurity(jParams, jAtomic);
            
            
            double Q1, Q2;
            jsx::array jDensityMatrix, jDensityMatrixPrime;
            
            read_density_matrix(jParams,
                                jAtomic,
                                jMeasurements,
                                Q1,
                                Q2,
                                jDensityMatrix,
                                jDensityMatrixPrime);
            
            
            /*
             std::cout << "Calculating impurity observables ...";
             if(jParams.is("OBS")) {
             throw std::runtime_error("static observables not yet implemented");
             
             jsx::value jObs; mpi::read(jParams("OBS").string(), jObs);
             
             for(auto& jEntry : jObs.object()) {
             double accObs = .0;
             for(std::size_t s = 0; s < jAtomic("Propagator").size(); ++s)
             accObs += linalg::trace(jsx::at<io::rmat>(jDensityMatrix.at(s)), jsx::at<io::rmat>(jEntry.second(s)));
             
             jScal[jEntry.first] = accObs;
             }
             }
             std::cout << "OK" << std::endl;
             */
        
            std::cout << "Calculating occupation ... ";
            
            {
                io::cmat occupation; occupation.resize(flavor_number, flavor_number);
                for(std::size_t i = 0; i < flavor_number; ++i) {
                    for(std::size_t j = 0; j < flavor_number; ++j) {
                        std::string entry = jParams("hybridisation")("matrix")(i)(j).string();
                        
                        if(!entry.empty())
                            occupation(i, j) = calculate_occupation(jDensityMatrix,
                                                                    jAtomic("operators dagger")(i),
                                                                    jAtomic("operators")(j));
                    }
                }
                auto occupation_entries = get_entries(occupation, jParams("hybridisation")("matrix"));
                
                for(auto const& entry : occupation_entries)
                    jObservables["occupation"][entry.first] = entry.second.real();
            }
            
            std::cout << "Ok" << std::endl;
            
            
            std::cout << "Calculating moments ... ";
            
            std::vector<io::cmat> sm(2);
            
            sm[0].resize(flavor_number, flavor_number);
            sm[1].resize(flavor_number, flavor_number);
            
            {
                std::vector<io::cmat> gm(3);
                
                gm[0].resize(flavor_number, flavor_number);
                gm[1].resize(flavor_number, flavor_number);
                gm[2].resize(flavor_number, flavor_number);
                
                for(std::size_t i = 0; i < flavor_number; ++i) {
                    for(std::size_t j = 0; j < flavor_number; ++j) {
                        std::string entry = jParams("hybridisation")("matrix")(i)(j).string();
                        
                        if(!entry.empty())
                            calculate_moments(jDensityMatrix,
                                              jDensityMatrixPrime,
                                              jAtomic("hloc")("eigen values"),
                                              jAtomic("operators")(i),
                                              jAtomic("operators dagger")(j),
                                              gm[1](i, j),
                                              gm[2](i, j));
                    }
                    
                    gm[0](i, i)  = 1.;
                    gm[1](i, i) += Q1;
                    gm[2](i, i) += Q2;
                }

                gm_entries.resize(3);
                
                gm_entries[0] = get_entries(gm[0], jParams("hybridisation")("matrix"));
                gm_entries[1] = get_entries(gm[1], jParams("hybridisation")("matrix"));
                gm_entries[2] = get_entries(gm[2], jParams("hybridisation")("matrix"));
                
                gm[1] = get_matrix(gm_entries[1], jParams("hybridisation")("matrix"));
                gm[2] = get_matrix(gm_entries[2], jParams("hybridisation")("matrix"));
                
                
                io::cmat gm1gm1; gm1gm1.resize(flavor_number, flavor_number);
                linalg::gemm('n', 'n', 1., gm[1], gm[1], .0, gm1gm1);
                
                for(std::size_t i = 0; i < flavor_number; ++i)
                    for(std::size_t j = 0; j < flavor_number; ++j) {
                        sm[0](i, j) += (i == j ? jParams("mu").real64() : .0) - jAtomic("hloc")("one body")(i)(j).real64() - gm[1](i, j);
                        sm[1](i, j) += gm[2](i, j) - gm1gm1(i, j);
                    }
                
                sm_entries.resize(2);
                
                sm_entries[0] = get_entries(sm[0], jParams("hybridisation")("matrix"));
                sm_entries[1] = get_entries(sm[1], jParams("hybridisation")("matrix"));
                
                sm[0] = get_matrix(sm_entries[0], jParams("hybridisation")("matrix")); //??
                sm[1] = get_matrix(sm_entries[1], jParams("hybridisation")("matrix")); //??
            }
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Fitting self-energy high frequency tail ... ";
            
            io::cmat alpha; alpha.resize(flavor_number, flavor_number);
            
            for(std::size_t i = 0; i < flavor_number; ++i)
                for(std::size_t j = 0; j < flavor_number; ++j)
                    if(jParams("hybridisation")("matrix")(i)(j).string() != "") {
                        std::size_t nFit  = std::max(self.size()/8, static_cast<std::size_t>(1));
                        std::complex<double> iomega_average = .0;
                        std::complex<double> self_average = .0;
                        
                        for(std::size_t n = self.size() - nFit;  n < self.size(); ++n) {
                            iomega_average += iomega(n);
                            self_average += self[n](i, j);
                        }
                        
                        alpha(i, j) = iomega_average/static_cast<double>(nFit) - sm[1](i, j)/(self_average/static_cast<double>(nFit) - sm[0](i, j));
                    }
            
            std::cout << "Ok" << std::endl;
            
            /*
            io::cmat alpha; alpha.resize(flavor_number, flavor_number);
            
            for(std::size_t i = 0; i < flavor_number; ++i)
                for(std::size_t j = 0; j < flavor_number; ++j)
                    if(jParams("hybridisation")("matrix")(i)(j).string() != "") {
                        std::size_t nFit  = std::max(self.size()/8, static_cast<std::size_t>(1));
                        
                        for(std::size_t n = self.size() - nFit;  n < self.size(); ++n)
                            alpha(i, j) -= sm[1](i, j)/(self[n](i, j) - sm[0](i, j)) - iomega(n);
                        
                        alpha(i, j) /= static_cast<double>(nFit);
                    }
            
            std::cout << "Ok" << std::endl;
            */
            
            std::cout << "Adding self-energy high frequency tail ... ";
            
            auto const omegaHF = iomega(self.size());
            
            for(std::size_t n = self.size(); n < hyb.size(); ++n) {
                io::cmat temp; temp.resize(flavor_number, flavor_number);
                
                for(std::size_t i = 0; i < flavor_number; ++i)
                    for(std::size_t j = 0; j < flavor_number; ++j)
                        temp(i, j) = sm[0](i, j) + sm[1](i, j)/(iomega(n) - std::complex<double>(alpha(i, j).real(), alpha(i, j).imag()*(omegaHF/iomega(n)).real()));
                
                self.push_back(temp);
            }
            
            std::cout << "OK" << std::endl;
            
            
            std::cout << "Adding green function high frequency tail ... ";
            
            for(std::size_t n = green.size(); n < hyb.size(); ++n) {
                io::cmat green_inv; green_inv.resize(flavor_number, flavor_number);
                
                for(std::size_t i = 0; i < flavor_number; ++i)
                    for(std::size_t j = 0; j < flavor_number; ++j)
                        green_inv(i, j) = (i == j ? iomega(n) + jParams("mu").real64() : .0) - jAtomic("hloc")("one body")(i)(j).real64() - hyb[n](i, j) - self[n](i, j);
                
                green.push_back(linalg::inv(green_inv));
            }
            
            std::cout << "OK" << std::endl;
            
        }
        
        
        {
            std::map<std::string, io::cvec> function_entries = get_function_entries(green, jParams("hybridisation")("matrix"));
            
            jsx::value jGreen;
            
            for(auto& function : function_entries) {
                jGreen[function.first]["function"] = std::move(function.second);
                
                if(gm_entries.size()) {
                    jGreen[function.first]["moments"] = jsx::array(3);
                    
                    jGreen[function.first]["moments"](0) =  gm_entries[0][function.first].real();
                    jGreen[function.first]["moments"](1) = -gm_entries[1][function.first].real();
                    jGreen[function.first]["moments"](2) =  gm_entries[2][function.first].real() + hm_entries[0][function.first].real();
                }
            }
            
            jObservables["green"] = std::move(jGreen);
        }
        
        
        {
            std::map<std::string, io::cvec> function_entries = get_function_entries(self, jParams("hybridisation")("matrix"));
            
            jsx::value jSelf;
            
            for(auto& function : function_entries) {
                jSelf[function.first]["function"] = std::move(function.second);
                
                if(sm_entries.size()) {
                    jSelf[function.first]["moments"] = jsx::array(2);
                    
                    jSelf[function.first]["moments"](0) = sm_entries[0][function.first].real();
                    jSelf[function.first]["moments"](1) = sm_entries[1][function.first].real();
                }
            }
            
            jObservables["self-energy"] =  std::move(jSelf);
        }
        
            
        mpi::write(jObservables, "observables.json");
        
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












