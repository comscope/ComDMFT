#ifndef LINALG_OPERATORS
#define LINALG_OPERATORS

#include "../io/Vector.h"
#include "../io/Matrix.h"
#include "LinAlg.h"


namespace linalg {
    
    jsx::value transpose(jsx::value& jOp) {
        jsx::value jOpDagg; jOpDagg = jsx::array_t(jOp.size(), jsx::object_t{{"target", jsx::null_t()}});
        
        int start_sector = 0;
        for(auto& jBlock : jOp.array()) {
            if(!jBlock("target").is<jsx::null_t>()) {
                int target_sector = jBlock("target").int64();
                
                jOpDagg[target_sector]["target"] = static_cast<jsx::int64_t>(start_sector);
                jOpDagg[target_sector]["matrix"] = jBlock("matrix");
                
                jsx::at<io::rmat>(jOpDagg[target_sector]["matrix"]).transpose();
            }
            ++start_sector;
        }
        
        return jOpDagg;
    }
    

    jsx::value diag_to_operator(jsx::value& jDiag) {
        jsx::value jOperator = jsx::array_t(jDiag.size());
        
        for(std::size_t s = 0; s < jDiag.size(); ++s) {
            std::size_t const dim = jsx::at<io::rvec>(jDiag(s)).size();
            
            io::rmat temp(dim, dim);
            for(std::size_t i = 0; i < dim; ++i)
                temp(i, i) = jsx::at<io::rvec>(jDiag(s))[i];
            
            jOperator(s)["matrix"] = std::move(temp);
            jOperator(s)["target"] = static_cast<jsx::int64_t>(s);
        }
        
        return jOperator;
    };

    
    void make_operator_complex(jsx::value& jOp) {
        for(auto& jBlock : jOp.array())
            if(!jBlock("target").is<jsx::null_t>()) {
                if(!jBlock.is<io::cmat>()) {
                    auto rmatrix = jsx::at<io::rmat>(jBlock("matrix"));
                    
                    io::cmat cmatrix(rmatrix.I(), rmatrix.J());
                    for(int i = 0; i < rmatrix.I(); ++i)
                        for(int j = 0; j < rmatrix.J(); ++j)
                            cmatrix(i, j) = rmatrix(i, j);
                    
                    jBlock("matrix") = cmatrix;
                }
            }
    };
    
    
    void mult(char transA, char transB, double alpha, jsx::value& jA, jsx::value& jB, double beta, jsx::value& jC)
    {
        if(transA != 'n' && transA != 't' && transA != 'c') throw std::runtime_error("im::mult: transA has invalid value");
        if(transB != 'n' && transB != 't' && transB != 'c') throw std::runtime_error("im::mult: transB has invalid value");
        if(jA.size() != jB.size()) throw std::runtime_error("im::mult: missmatch in sector number of A and B");
        
        if(beta == .0) jC = jsx::array_t(jA.size(), jsx::object_t{{"target", jsx::null_t()}});
        
        jsx::value jMapA = jsx::array_t(jA.size(), jsx::null_t());
        jsx::value jMapB = jsx::array_t(jA.size(), jsx::null_t());
        
        for(std::size_t s = 0; s < jA.size(); ++s) {
            if(!jA(s)("target").is<jsx::null_t>())
                transA == 'n' ? jMapA(s) = jA(s)("target").int64() : jMapA(jA(s)("target").int64()) = jsx::int64_t(s);
            
            if(!jB(s)("target").is<jsx::null_t>())
                transB == 'n' ? jMapB(s) = jB(s)("target").int64() : jMapB(jB(s)("target").int64()) = jsx::int64_t(s);
        }
        
        std::vector<bool> guard(jA.size(), false);
        
        for(std::size_t s = 0; s < jMapB.size(); ++s) {
            if(!jMapB(s).is<jsx::null_t>()) {
                std::size_t const m = jMapB(s).int64();
                
                if(!jMapA(m).is<jsx::null_t>()) {
                    std::size_t const t = jMapA(m).int64();
                    
                    if(guard.at(t))
                        throw std::runtime_error("im::mult: A or B has invalid block structure");
                    
                    int const secA = transA == 'n' ? m : t; int const secB = transB == 'n' ? s : m;
                    
                    if(jA(secA)("matrix").is<io::rmat>() && jB(secB)("matrix").is<io::rmat>()) {
                        io::rmat const& matA = jsx::at<io::rmat>(jA(secA)("matrix"));
                        io::rmat const& matB = jsx::at<io::rmat>(jB(secB)("matrix"));
                        
                        if(jC(s)("target").is<jsx::null_t>()) {
                            jC(s)("target") = jsx::int64_t(t);
                            jC(s)["matrix"] = io::rmat(transA == 'n' ? matA.I() : matA.J(), transB == 'n' ? matB.J() : matB.I());
                        }
                        
                        if(jC(s)("target").int64() != jsx::int64_t(t))
                            throw std::runtime_error("im::mult: block structure of C not compatible with A and B");
                        
                        linalg::gemm(transA == 'c' ? 't' : transA, transB == 'c' ? 't' : transB, alpha, matA, matB, beta, jsx::at<io::rmat>(jC(s)("matrix")));
                    } else if(jA(secA)("matrix").is<io::cmat>() && jB(secB)("matrix").is<io::cmat>()) {
                        io::cmat const& matA = jsx::at<io::cmat>(jA(secA)("matrix"));
                        io::cmat const& matB = jsx::at<io::cmat>(jB(secB)("matrix"));
                        
                        if(jC(s)("target").is<jsx::null_t>()) {
                            jC(s)("target") = jsx::int64_t(t);
                            jC(s)["matrix"] = io::cmat(transA == 'n' ? matA.I() : matA.J(), transB == 'n' ? matB.J() : matB.I());
                        }
                        
                        if(jC(s)("target").int64() != jsx::int64_t(t))
                            throw std::runtime_error("im::mult: block structure of C not compatible with A and B");
                        
                        linalg::gemm(transA, transB, alpha, matA, matB, beta, jsx::at<io::cmat>(jC(s)("matrix")));
                    } else
                        throw std::runtime_error("linalg::mult: unknow matrix format or matrix type combination");
                    
                    guard.at(t) = true;
                }
            }
        }
    };
    
    
    std::complex<double> trace(jsx::value& jA, jsx::value& jB)
    {
        if(jA.size() != jB.size()) throw std::runtime_error("im::mult: missmatch in sector number of A and B");
        
        std::complex<double> result = .0;
        
        for(std::size_t s = 0; s < jA.size(); ++s)
            if(!jA(s)("target").is<jsx::null_t>()) {
                std::size_t const m = jA(s)("target").int64();
                
                if(!jB(m)("target").is<jsx::null_t>()) {
                    std::size_t const t = jB(m)("target").int64();
                    
                    if(s != t)
                        throw std::runtime_error("im::trace: product is not block-diagonal");
                    
                    if(jA(s)("matrix").is<io::rmat>() && jB(s)("matrix").is<io::rmat>())
                        
                        result += linalg::trace(jsx::at<io::rmat>(jA(s)("matrix")), jsx::at<io::rmat>(jB(s)("matrix")));
                    
                    else if(jA(s)("matrix").is<io::cmat>() && jB(s)("matrix").is<io::cmat>())
                        
                        result += linalg::trace(jsx::at<io::cmat>(jA(s)("matrix")), jsx::at<io::cmat>(jB(s)("matrix")));
                    
                    else
                        
                        throw std::runtime_error("linalg::trace: unknow matrix format or matrix type combination");
                }
            }
        
        return result;
    };
}

#endif










