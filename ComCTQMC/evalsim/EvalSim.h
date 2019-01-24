#ifndef EVALSIM
#define EVALSIM

#include "../include/Measurements.h"
#include "../include/LinAlg.h"
#include "../include/IO.h"
#include "../include/impurity/GenerateAtomic.h"

inline void setup_impurity(jsx::value const& jParams, jsx::value& jAtomic) {
    std::cout << "Reading in eigen values ...";
    
    auto& jEigenValues = jAtomic("hloc")("eigen values");
    auto& jN = jAtomic("hloc")("quantum numbers")("N");
    for(std::size_t sector = 0; sector < jEigenValues.size(); ++sector)
        for(auto& energy : jsx::at<io::rvec>(jEigenValues(sector)))
            energy += -jParams("mu").real64()*jsx::at<io::rvec>(jN)[sector];
    
    std::cout << "Ok" << std::endl;
    
    
    std::cout << "Reading in operators ...";
    
    jsx::array jOperatorsDagger;
    for(auto& jOp : jAtomic("operators").array()) {
        jsx::value jOpDagger; jOpDagger = jsx::array(jOp.size(), jsx::object{{"target", jsx::null()}});
        
        int start_sector = 0;
        for(auto& jBloc : jOp.array()) {
            if(jBloc("target").type() != jsx::type::null) {
                int target_sector = jBloc("target").int64();
                
                jOpDagger[target_sector]["target"] = static_cast<std::int64_t>(start_sector);
                jOpDagger[target_sector]["matrix"] = jBloc("matrix");
                
                jsx::at<io::rmat>(jOpDagger[target_sector]["matrix"]).transpose();
            }
            ++start_sector;
        }
        jOperatorsDagger.push_back(std::move(jOpDagger));
    }
    jAtomic["operators dagger"] = std::move(jOperatorsDagger);
    
    std::cout << "Ok" << std::endl;
}


inline void read_density_matrix(jsx::value const& jParams,
                         jsx::value& jAtomic,
                         jsx::value const& jMeasurements,
                         double& Q1,
                         double& Q2,
                         jsx::array& jDensityMatrix,
                         jsx::array& jDensityMatrixPrime) {
    
    std::cout << "Reading in density matrix ...";
    
    Q1 = Q2 = .0;
    
    double UrOmega0 = .0;
    
    if(jParams.is("dyn")) {
        jsx::value jDynF0; mpi::read(jParams("dyn").string(), jDynF0);
        
        UrOmega0 = jDynF0(0).real64();
        
        double UrTau0 = .0;
        for(std::size_t n = 0; n < jDynF0.size(); ++n)
            UrTau0 += 2./jParams("beta").real64()*jDynF0(n).real64();
        
        Q1 = -UrOmega0*jsx::at<io::rvec>(jMeasurements("Scal")("N"))[0];
        Q2 = -UrTau0 + jsx::at<io::rvec>(jMeasurements("DynDyn"))[0];
    }
    
    for(std::size_t s = 0; s < jAtomic("hloc")("eigen values").size(); ++s) {
        int const dim = jsx::at<io::rvec>(jAtomic("hloc")("eigen values")(s)).size(); /// at not necessary ??
        
        auto const& source = jsx::at<io::rvec>(jMeasurements("DensityMatrix")(s));
        
        io::rmat target; target.resize(dim, dim);
        for(int i = 0; i < dim; ++i)
            for(int j = 0; j < dim; ++j)
                target(i, j) = (source.at(i*dim + j) + source.at(j*dim + i))/2.;
        
        io::rmat targetPrime; targetPrime.resize(dim, dim);
        if(jParams.is("dyn")) {
            int const filling = jsx::at<io::rvec>(jAtomic("hloc")("quantum numbers")("N"))[s];
            
            auto const& sourceDyn = jsx::at<io::rvec>(jMeasurements("DensityMatrixDyn")(s));
            
            for(int i = 0; i < dim; ++i)
                for(int j = 0; j < dim; ++j)
                    targetPrime(i, j) = UrOmega0*filling*target(i, j) + (sourceDyn.at(i*dim + j) + sourceDyn.at(j*dim + i))/2.;
            
            Q2 += 2.*UrOmega0*filling*linalg::trace(targetPrime) - UrOmega0*UrOmega0*filling*filling*linalg::trace(target);
        }
        
        jDensityMatrix.push_back(std::move(target));
        jDensityMatrixPrime.push_back(std::move(targetPrime));
    }
    
    std::cout << "Ok" << std::endl;
}


inline double calculate_occupation(jsx::array& jDensityMatrix, jsx::value& jOpDagg, jsx::value& jOp) {
    double occupation = .0;
    
    for(std::size_t s = 0; s < jDensityMatrix.size(); ++s) {
        if(jOp(s)("target").type() != jsx::type::null) {
            std::size_t const m = jOp(s)("target").int64();
            
            if(jOpDagg(m)("target").type() != jsx::type::null) {
                std::size_t const t = jOpDagg(m)("target").int64();
                
                if(t == s) {
                    auto const& opBloc = jsx::at<io::rmat>(jOp(s)("matrix"));
                    auto const& opDaggBloc = jsx::at<io::rmat>(jOpDagg(m)("matrix"));
                    
                    io::rmat temp; temp.resize(opDaggBloc.I(), opBloc.J());
                    linalg::gemm('n', 'n', 1., opDaggBloc, opBloc, .0, temp);
                    occupation += linalg::trace(jsx::at<io::rmat>(jDensityMatrix.at(s)), temp);
                } else
                    throw std::runtime_error("occupation can not be calculated from reduced density matrix");
            }
        }
    }
    
    return occupation;
};


inline void calculate_moments(jsx::array& jDensityMatrix,
                       jsx::array& jDensityMatrixPrime,
                       jsx::value& jEigenValues,
                       jsx::value& jOp,
                       jsx::value& jOpDagg,
                       std::complex<double>& g1, std::complex<double>& g2)
{
    double temp1 = .0; double temp2 = .0;
    
    for(std::size_t s = 0; s < jEigenValues.size(); ++s) {
        //~ op opDagg
        if(jOpDagg(s)("target").type() != jsx::type::null) {
            std::size_t const m = jOpDagg(s)("target").int64();
            
            if(jOp(m)("target").type() != jsx::type::null) {
                std::size_t const t = jOp(m)("target").int64();
                
                if(t == s) {
                    auto const& opDaggBloc = jsx::at<io::rmat>(jOpDagg(s)("matrix"));
                    auto const& opBloc = jsx::at<io::rmat>(jOp(m)("matrix"));
                    
                    auto const& sHlocBloc = jsx::at<io::rvec>(jEigenValues(s));
                    auto const& mHlocBloc = jsx::at<io::rvec>(jEigenValues(m));
                    auto const& tHlocBloc = jsx::at<io::rvec>(jEigenValues(t));
                    
                    //C=[H,op]
                    io::rmat C; C.resize(opBloc.I(), opBloc.J());
                    
                    for(int a = 0; a < opBloc.I(); ++a)
                        for(int b = 0; b < opBloc.J(); ++b)
                            C(a, b) = (tHlocBloc[a] - mHlocBloc[b])*opBloc(a, b);
                    
                    //CDagg=[opDagg,H]
                    io::rmat CDagg; CDagg.resize(opDaggBloc.I(), opDaggBloc.J());
                    
                    for(int a = 0; a < opDaggBloc.I(); ++a)
                        for(int b = 0; b < opDaggBloc.J(); ++b)
                            CDagg(a, b) = opDaggBloc(a, b)*(sHlocBloc[b] - mHlocBloc[a]);
                    
                    io::rmat temp; temp.resize(opBloc.I(), opDaggBloc.J());
                    
                    linalg::gemm('n', 'n', 1., C, opDaggBloc, .0, temp);
                    temp1 += linalg::trace(jsx::at<io::rmat>(jDensityMatrix.at(t)), temp);
                    
                    linalg::gemm('n', 'n', 1., C, CDagg, .0, temp);
                    temp2 += linalg::trace(jsx::at<io::rmat>(jDensityMatrix.at(t)), temp);
                    
                    linalg::gemm('n', 'n', 1., opBloc, CDagg, .0, temp);
                    temp2 -= linalg::trace(jsx::at<io::rmat>(jDensityMatrixPrime.at(t)), temp);
                    
                    linalg::gemm('n', 'n', 1., C, opDaggBloc, .0, temp);
                    temp2 -= linalg::trace(jsx::at<io::rmat>(jDensityMatrixPrime.at(t)), temp);
                } else
                    throw std::runtime_error("Moment can not be calculated from reduced density matrix");                                }
        }
        
        //~ opDagg op
        if(jOp(s)("target").type() != jsx::type::null) {
            std::size_t const m = jOp(s)("target").int64();
            
            if(jOpDagg(m)("target").type() != jsx::type::null) {
                std::size_t const t = jOpDagg(m)("target").int64();
                
                if(t == s) {
                    auto const& opBloc = jsx::at<io::rmat>(jOp(s)("matrix"));
                    auto const& opDaggBloc = jsx::at<io::rmat>(jOpDagg(m)("matrix"));
                    
                    auto const& sHlocBloc = jsx::at<io::rvec>(jEigenValues(s));
                    auto const& mHlocBloc = jsx::at<io::rvec>(jEigenValues(m));
                    auto const& tHlocBloc = jsx::at<io::rvec>(jEigenValues(t));
                    
                    //C=[H,op]
                    io::rmat C; C.resize(opBloc.I(), opBloc.J());
                    
                    for(int a = 0; a < opBloc.I(); ++a)
                        for(int b = 0; b < opBloc.J(); ++b)
                            C(a, b) = (mHlocBloc[a] - sHlocBloc[b])*opBloc(a, b);
                    
                    //CDagg=[opDagg,H]
                    io::rmat CDagg; CDagg.resize(opDaggBloc.I(), opDaggBloc.J());
                    
                    for(int a = 0; a < opDaggBloc.I(); ++a)
                        for(int b = 0; b < opDaggBloc.J(); ++b)
                            CDagg(a, b) = opDaggBloc(a, b)*(mHlocBloc[b] - tHlocBloc[a]);
                    
                    io::rmat temp; temp.resize(opDaggBloc.I(), opBloc.J());
                    
                    linalg::gemm('n', 'n', 1., opDaggBloc, C, .0, temp);
                    temp1 += linalg::trace(jsx::at<io::rmat>(jDensityMatrix.at(t)), temp);
                    
                    linalg::gemm('n', 'n', 1., CDagg, C, .0, temp);
                    temp2 += linalg::trace(jsx::at<io::rmat>(jDensityMatrix.at(t)), temp);
                    
                    linalg::gemm('n', 'n', 1., CDagg, opBloc, .0, temp);
                    temp2 -= linalg::trace(jsx::at<io::rmat>(jDensityMatrixPrime.at(t)), temp);
                    
                    linalg::gemm('n', 'n', 1., opDaggBloc, C, .0, temp);
                    temp2 -= linalg::trace(jsx::at<io::rmat>(jDensityMatrixPrime.at(t)), temp);
                } else
                    throw std::runtime_error("Moment can not be calculated from reduced density matrix");
            }
        }
    }
    
    g1 = temp1; g2 = temp2;
}

io::cmat get_matrix(std::map<std::string, std::complex<double>> const& entries, jsx::value const& jMatrix) {
    io::cmat matrix; matrix.resize(jMatrix.size(), jMatrix.size());
    
    for(std::size_t i = 0; i < jMatrix.size(); ++i)
        for(std::size_t j = 0; j < jMatrix.size(); ++j) {
            auto const entry = jMatrix(i)(j).string();
            
            if(entry != "") matrix(i, j) = entries.at(entry);
        }
    
    return matrix;
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
};


inline std::vector<io::cmat> get_function_matrix(std::map<std::string, io::cvec> const& functions, jsx::value const& jMatrix)
{
    std::size_t size = std::numeric_limits<std::size_t>::max();
    
    for(auto const& function : functions)
        size = std::min(size, function.second.size());
    
    std::vector<io::cmat> function_matrix(size);
    for(std::size_t n = 0; n < size; ++n) {
        auto& matrix = function_matrix[n];
        
        matrix.resize(jMatrix.size(), jMatrix.size());
        for(std::size_t i = 0; i < jMatrix.size(); ++i)
            for(std::size_t j = 0; j < jMatrix.size(); ++j) {
                auto const& entry = jMatrix(i)(j).string();
                
                if(entry != "")
                    matrix(i, j) = functions.at(entry)[n];
            }
    }
    
    return function_matrix;
}

inline std::map<std::string, io::cvec> get_function_entries(std::vector<io::cmat> const& function_matrix, jsx::value const& jMatrix)
{
    std::map<std::string, io::cvec> function_entries;
    
    for(auto const& matrix : function_matrix) {
        std::map<std::string, std::complex<double>> entries = get_entries(matrix, jMatrix);
        for(auto const& entry : entries)
            function_entries[entry.first].push_back(entry.second);
    }
    
    return function_entries;
}


struct iOmega {
    iOmega() = delete;
    iOmega(double beta) : beta_(beta) {};
    std::complex<double> operator()(int n) {
        return {.0, (2*n + 1)*M_PI/beta_};
    };
private:
    double const beta_;
};

#endif //EVALSIM










