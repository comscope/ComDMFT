#ifndef GENERATEATOMIC
#define GENERATEATOMIC

#include <cmath>
#include <stdexcept>
#include <iostream>
#include <vector>
#include <bitset>
#include <cassert>
#include <iomanip>
#include <set>
#include <numeric>
#include <algorithm>
#include <fstream>
#include <string>

#include "../JsonX.h"
#include "../LinAlg.h"
#include "../IO.h"

namespace Ga {
    //kack isch veralted das muess aapassed werde ....  c++11 !!!!!!
    
    struct FlavorState : public std::bitset<8*sizeof(unsigned long long)> {
        typedef unsigned long long State;
        
        FlavorState() : sign_(1) {};
        explicit FlavorState(State const& state) : std::bitset<8*sizeof(State)>(state), sign_(1) {};
        FlavorState(FlavorState const& state) : std::bitset<8*sizeof(State)>(state), sign_(state.sign()) {};
        
        State state() const { return this->to_ullong();};
        int sign() const { return sign_;};
        int& sign() { return sign_;};
    private:
        int sign_;
    };
    
    
    FlavorState psi(int flavor, FlavorState const& state) {
        FlavorState temp = state;
        if(state.test(flavor)) {
            temp[flavor] = 0; if((state >> (flavor + 1)).count()%2) temp.sign() *= -1;
        } else
            temp.sign() = 0;
        return temp;
    };
    
    FlavorState psiDagg(int flavor, FlavorState const& state) {
        FlavorState temp = state;
        if(state.test(flavor))
            temp.sign() = 0;
        else {
            temp[flavor] = 1; if((state >> (flavor + 1)).count()%2) temp.sign() *= -1;
        }
        return temp;
    };
    
    typedef FlavorState::State State;
    typedef std::vector<State> States;
    typedef std::vector<States> BlocStates;
    
    
    struct OperatorBloc {
        int& target() { return target_;};
        int target() const { return target_;};
        
        io::rmat& matrix() { return matrix_;};
        io::rmat const& matrix() const { return matrix_;};
    private:
        int target_ = -1; ////!!!!!!!!!!!!!!!!!!!!!!!!
        io::rmat matrix_;
    };
    
    
    
    
    
    template<typename Hloc>
    void find_blocks(Hloc const& hloc, BlocStates& blocStates)
    {
        std::set<State> states;
        for(State u = 0; u < static_cast<State>(1 << hloc.N()); ++u) states.insert(u);
        
        while(states.size()) {
            blocStates.push_back(States());
            
            // implement cluster labeling algorithm as below
            std::set<State> temp; temp.insert(*states.begin());
            while(temp.size()) {
                FlavorState const startState(*temp.begin()); temp.erase(temp.begin());
                blocStates.back().push_back(startState.state());
                
                for(int fDagg = 0; fDagg < hloc.N(); ++fDagg)
                    for(int f = 0; f < hloc.N(); ++f) {
                        FlavorState const targetState = psiDagg(fDagg, psi(f, startState));
                        
                        if((targetState.sign() != 0) && (hloc.t(fDagg, f) != .0))
                            if(std::find(blocStates.back().begin(), blocStates.back().end(), targetState.state()) == blocStates.back().end())
                                temp.insert(targetState.state());
                    }
                
                for(int f1Dagg = 0; f1Dagg < hloc.N(); ++f1Dagg)
                    for(int f2Dagg = 0; f2Dagg < hloc.N(); ++f2Dagg)
                        for(int f1 = 0; f1 < hloc.N(); ++f1)
                            for(int f2 = 0; f2 < hloc.N(); ++f2) {
                                FlavorState const targetState = psiDagg(f1Dagg, psiDagg(f2Dagg, psi(f1, psi(f2, startState))));
                                
                                if((targetState.sign() != 0) && (hloc.V(f1Dagg, f2Dagg, f1, f2) != .0))
                                    if(std::find(blocStates.back().begin(), blocStates.back().end(), targetState.state()) == blocStates.back().end())
                                        temp.insert(targetState.state());
                            }
                
                states.erase(startState.state());
            }
        }
        
        std::vector<std::size_t> index(blocStates.size()); std::iota(index.begin(), index.end(), 0);
        
        for(BlocStates::const_iterator itBlocStatesJ = blocStates.begin(); itBlocStatesJ != blocStates.end(); ++itBlocStatesJ)
            for(int f = 0; f < hloc.N(); ++f)
                for(int type = 0; type < 2; ++type) {
                    int targetSector = -1;
                    for(States::const_iterator itStatesJ = itBlocStatesJ->begin(); itStatesJ != itBlocStatesJ->end(); ++itStatesJ) {
                        FlavorState const stateJ(*itStatesJ);
                        FlavorState const stateI = type ? psiDagg(f, stateJ) : psi(f, stateJ);;
                        if(stateI.sign() != 0) {
                            int targetSectorTemp = 0;
                            BlocStates::const_iterator itBlocStatesI = blocStates.begin();
                            for(; std::find(itBlocStatesI->begin(), itBlocStatesI->end(), stateI.state()) == itBlocStatesI->end(); ++itBlocStatesI, ++targetSectorTemp);
                            if(targetSector < 0)
                                targetSector = targetSectorTemp;
                            else if(targetSector != targetSectorTemp) {
                                std::size_t i = index[targetSector];
                                while(i != index[i]) i = index[i];
                                
                                std::size_t iTemp = index[targetSectorTemp];
                                while(iTemp != index[iTemp]) iTemp = index[iTemp];
                                
                                index[iTemp] < index[i] ? index[i] = index[iTemp] : index[iTemp] = index[i];
                            }
                        }
                    }
                }
        
        BlocStates blocStatesTemp(blocStates.size());
        for(std::size_t sector = 0; sector < blocStates.size(); ++sector) {
            std::size_t i = index[sector];
            while(i != index[i]) i = index[i];
            blocStatesTemp[i].insert(blocStatesTemp[i].end(), blocStates[sector].begin(), blocStates[sector].end());
        };
        
        blocStates.clear();
        for(std::size_t sector = 0; sector < blocStatesTemp.size(); ++sector)
            if(blocStatesTemp[sector].size()) blocStates.push_back(blocStatesTemp[sector]);
    };
    
    
    template<typename Hloc>
    void construct_hamiltonian(Hloc const& hloc,
                               BlocStates const& blocStates,
                               std::vector<io::rmat>& hamiltonian)
    {
        hamiltonian.resize(blocStates.size());
        
        int startSector = 0;
        for(BlocStates::const_iterator itBlocStatesJ = blocStates.begin(); itBlocStatesJ != blocStates.end(); ++itBlocStatesJ, ++startSector) {
            hamiltonian[startSector].resize(itBlocStatesJ->size(), itBlocStatesJ->size());
            
            int indexJ = 0;
            for(States::const_iterator itStatesJ = itBlocStatesJ->begin(); itStatesJ != itBlocStatesJ->end(); ++itStatesJ, ++indexJ) {
                FlavorState const stateJ(*itStatesJ);
                
                for(int fDagg = 0; fDagg < hloc.N(); ++fDagg)
                    for(int f = 0; f < hloc.N(); ++f) {
                        FlavorState const stateI = psiDagg(fDagg, psi(f, stateJ));
                        
                        if(stateI.sign() != 0) {
                            double const t = hloc.t(fDagg, f);
                            
                            if(t != .0) {
                                auto itIndexI = std::find(itBlocStatesJ->begin(), itBlocStatesJ->end(), stateI.state());
                                
                                if(itIndexI == itBlocStatesJ->end())
                                    throw std::runtime_error("Something is wrong with the partitioning of the states");
                                
                                hamiltonian[startSector](itIndexI - itBlocStatesJ->begin(), indexJ) += t*stateI.sign();
                            }
                        }
                    }
                
                for(int f1Dagg = 0; f1Dagg < hloc.N(); ++f1Dagg)
                    for(int f2Dagg = 0; f2Dagg < hloc.N(); ++f2Dagg)
                        for(int f1 = 0; f1 < hloc.N(); ++f1)
                            for(int f2 = 0; f2 < hloc.N(); ++f2) {
                                FlavorState const stateI = psiDagg(f1Dagg, psiDagg(f2Dagg, psi(f1, psi(f2, stateJ))));
                                
                                if(stateI.sign() != 0) {
                                    double const V = hloc.V(f1Dagg, f2Dagg, f1, f2);
                                    
                                    if(V != .0) {
                                        auto itIndexI = std::find(itBlocStatesJ->begin(), itBlocStatesJ->end(), stateI.state());
                                        
                                        if(itIndexI == itBlocStatesJ->end())
                                            throw std::runtime_error("Something is wrong with the partitioning of the states");
                                        
                                        hamiltonian[startSector](itIndexI - itBlocStatesJ->begin(), indexJ) += V*stateI.sign();
                                    }
                                }
                            }
            }
        }
    };
    
    
    
    void diagonalise(std::vector<io::rmat>& hamiltonian,
                     std::vector<io::rvec>& eigen_values)
    {
        eigen_values.resize(hamiltonian.size());
        
        for(unsigned int sector = 0; sector < hamiltonian.size(); ++sector) {
            eigen_values[sector].resize(hamiltonian[sector].I());
            linalg::syev('V', 'U', hamiltonian[sector], eigen_values[sector]);
        }
    };
    
    
    
    void construct_operators(int N,
                             BlocStates const& blocStates,
                             std::vector<std::vector<OperatorBloc>>& operators)
    {
        operators.resize(N, std::vector<OperatorBloc>(blocStates.size())); // c.f. definition of OperatorBloc:  target = -1
        
        int startSector = 0;
        for(BlocStates::const_iterator itBlocStatesJ = blocStates.begin(); itBlocStatesJ != blocStates.end(); ++itBlocStatesJ, ++startSector) {
            int indexJ = 0;
            for(States::const_iterator itStatesJ = itBlocStatesJ->begin(); itStatesJ != itBlocStatesJ->end(); ++itStatesJ, ++indexJ) {
                FlavorState const stateJ(*itStatesJ);
                for(int f = 0; f < N; ++f)  {
                    FlavorState const stateI = psi(f, stateJ);
                    if(stateI.sign() != 0) {
                        int targetSector = 0;
                        States::const_iterator itStatesI;
                        BlocStates::const_iterator itBlocStatesI = blocStates.begin();
                        for(; (itStatesI = std::find(itBlocStatesI->begin(), itBlocStatesI->end(), stateI.state())) == itBlocStatesI->end(); ++itBlocStatesI, ++targetSector);
                        int indexI = itStatesI - itBlocStatesI->begin();
                        
                        if(operators[f][startSector].target() == -1) {
                            operators[f][startSector].target() = targetSector;
                            operators[f][startSector].matrix().resize(itBlocStatesI->size(), itBlocStatesJ->size());
                        } else if(operators[f][startSector].target() != targetSector)
                            throw std::runtime_error("Fatal error: Operator maps one sector to two or more sectors.");
                        
                        operators[f][startSector].matrix()(indexI, indexJ) += stateI.sign();
                    }
                }
            }
        }
    };
    
    
    
    void transform_operators(std::vector<io::rmat> const& transformation,
                             std::vector<std::vector<OperatorBloc>>& operators)
    {
        for(int f = 0; f < static_cast<int>(operators.size()); ++f)
            for(int sector = 0; sector < static_cast<int>(transformation.size()); ++sector) {
                if(operators[f][sector].target() != -1) {
                    io::rmat buffer; buffer.resize(operators[f][sector].matrix().I(), operators[f][sector].matrix().J());
                    linalg::gemm('n', 'n', 1., operators[f][sector].matrix(), transformation[sector], .0, buffer);
                    linalg::gemm('t', 'n', 1., transformation[operators[f][sector].target()], buffer, .0, operators[f][sector].matrix());
                }
            }
    };
    
    
    
    void calculate_QN(BlocStates const& blocStates, std::vector<double> const& qn, io::rvec& QN)
    {
        for(BlocStates::const_iterator itBloc = blocStates.begin(); itBloc != blocStates.end(); ++itBloc)  {
            double value = .0;
            
            for(std::size_t f = 0; f < qn.size(); ++f)
                if(FlavorState(*itBloc->begin()).test(f)) value += qn[f];
            
            for(States::const_iterator itStates = itBloc->begin(); itStates != itBloc->end(); ++itStates) {
                double temp = .0;
                
                for(std::size_t f = 0; f < qn.size(); ++f)
                    if(FlavorState(*itStates).test(f)) temp += qn[f];
                
                if(value != temp)
                    throw std::runtime_error("Problem with quantum numbers.");
            }
            
            QN.push_back(value);
        };
    };
    
    
    template<typename Obs>
    void construct_observable(Obs const& obs,
                              BlocStates const& blocStates,
                              std::vector<io::rmat>& observable)
    {
        observable.resize(blocStates.size());
        
        int startSector = 0;
        for(BlocStates::const_iterator itBlocStatesJ = blocStates.begin(); itBlocStatesJ != blocStates.end(); ++itBlocStatesJ, ++startSector) {
            observable[startSector].resize(itBlocStatesJ->size(), itBlocStatesJ->size());
            
            int indexJ = 0;
            for(States::const_iterator itStatesJ = itBlocStatesJ->begin(); itStatesJ != itBlocStatesJ->end(); ++itStatesJ, ++indexJ) {
                FlavorState const stateJ(*itStatesJ);
                
                for(int fDagg = 0; fDagg < obs.N(); ++fDagg)
                    for(int f = 0; f < obs.N(); ++f) {
                        FlavorState const stateI = psiDagg(fDagg, psi(f, stateJ));
                        
                        if(stateI.sign() != 0) {
                            double const t = obs.t(fDagg, f);
                            
                            if(t != .0) {
                                auto itIndexI = std::find(itBlocStatesJ->begin(), itBlocStatesJ->end(), stateI.state());
                                
                                if(itIndexI == itBlocStatesJ->end())
                                    throw std::runtime_error("Matrix observable does not respect the symmetries");
                                
                                observable[startSector](itIndexI - itBlocStatesJ->begin(), indexJ) += t*stateI.sign();
                            }
                        }
                    }
                
                for(int f1Dagg = 0; f1Dagg < obs.N(); ++f1Dagg)
                    for(int f1 = 0; f1 < obs.N(); ++f1)
                        for(int f2Dagg = 0; f2Dagg < obs.N(); ++f2Dagg)
                            for(int f2 = 0; f2 < obs.N(); ++f2) {
                                FlavorState const stateI = psiDagg(f1Dagg, psi(f1, psiDagg(f2Dagg, psi(f2, stateJ))));
                                
                                if(stateI.sign() != 0) {
                                    double const V = obs.V(f1Dagg, f2Dagg, f1, f2);
                                    
                                    if(V != .0) {
                                        auto itIndexI = std::find(itBlocStatesJ->begin(), itBlocStatesJ->end(), stateI.state());
                                        
                                        if(itIndexI == itBlocStatesJ->end())
                                            throw std::runtime_error("Matrix observable does not respect the symmetries");
                                        
                                        observable[startSector](itIndexI - itBlocStatesJ->begin(), indexJ) += V*stateI.sign();
                                    }
                                }
                            }
            }
        }
    };
    
    
    void transform_observable(std::vector<io::rmat> const& transformation,
                              std::vector<io::rmat>& observable)
    {
        for(int sector = 0; sector < static_cast<int>(transformation.size()); ++sector) {
            io::rmat buffer; buffer.resize(transformation[sector].I(), transformation[sector].J());
            linalg::gemm('n', 'n', 1., observable[sector], transformation[sector], .0, buffer);
            linalg::gemm('t', 'n', 1., transformation[sector], buffer, .0, observable[sector]);
        }
    };
    
    
    struct GenerateAtomic {
        GenerateAtomic() = delete;
        
        template<typename Hloc>
        GenerateAtomic(Hloc const& hloc, jsx::value& jHloc, bool b64 = true) : flavor_number_(hloc.N()) {
            
            std::cout << "Number of invariant subspaces: ";
            find_blocks(hloc, blocStates_);
            std::cout << blocStates_.size() << std::endl;
            
            std::vector<io::rmat> hamiltonian;
            construct_hamiltonian(hloc, blocStates_, hamiltonian);
            
            int maxDim = 0;
            for(auto& block : hamiltonian)
                maxDim = std::max(block.I(), maxDim);
            
            std::cout << "Dimension of the biggest subspace: " << maxDim << std::endl;
            
            std::vector<io::rvec> eigen_values;
            diagonalise(hamiltonian, eigen_values);
            transformation_ = std::move(hamiltonian);
            
            std::map<std::string, io::rvec> QNS;
            for(auto const& qn : hloc.qns())
                calculate_QN(blocStates_, qn.second, QNS[qn.first]);
            
            jHloc["flavor number"] = static_cast<std::int64_t>(flavor_number_);
            
            jsx::array jBlocStates;
            for(auto const& states : blocStates_) {
                jsx::array jStates;
                for(auto const& state : states)
                    jStates.push_back(static_cast<std::int64_t>(state));
                jBlocStates.push_back(std::move(jStates));
            }
            jHloc["block states"] = std::move(jBlocStates);
            
            jsx::array jTransformation;
            for(auto const& bloc : transformation_) {
                bloc.b64() = b64;
                jTransformation.push_back(bloc);
            }
            jHloc["transformation"] = std::move(jTransformation);
            
            jsx::array jEigenValues;
            for(auto const& bloc : eigen_values) {
                bloc.b64() = b64;
                jEigenValues.push_back(bloc);
            }
            jHloc["eigen values"] = std::move(jEigenValues);
            
            for(auto& QN : QNS) {
                QN.second.b64() = b64;
                jHloc["quantum numbers"][QN.first] = std::move(QN.second);
            }
            
        };
        
        
        GenerateAtomic(jsx::value& jHloc) : flavor_number_(jHloc("flavor number").int64()) {
            if(jHloc("block states").size() != jHloc("transformation").size())
                throw std::runtime_error("GenerateAtomic: sector number differs");
            
            for(std::size_t sector = 0; sector < jHloc("block states").size(); ++sector) {
                auto const& block = jsx::at<io::rmat>(jHloc("transformation")(sector));
                
                if(jHloc("block states")(sector).size() != static_cast<std::size_t>(block.I()) || jHloc("block states")(sector).size() != static_cast<std::size_t>(block.J()))
                    throw std::runtime_error("GenerateAtomic: transformation and eigen values not compatible");
            }
            
            for(auto const& jStates : jHloc("block states").array()) {
                States states;
                for(auto const& jState : jStates.array())
                    states.push_back(jState.int64());
                blocStates_.push_back(std::move(states));
            }
            
            for(auto& jBloc : jHloc("transformation").array())
                transformation_.push_back(jsx::at<io::rmat>(jBloc));
        };
        
        
        GenerateAtomic(GenerateAtomic const&) = delete;
        GenerateAtomic(GenerateAtomic&&) = delete;
        GenerateAtomic& operator=(GenerateAtomic const&) = delete;
        GenerateAtomic& operator=(GenerateAtomic&&) = delete;
        ~GenerateAtomic() = default;
        
        void write_operators(jsx::value& jOperators, int b64 = true) const {
            
            std::vector<std::vector<OperatorBloc> > operators;
            construct_operators(flavor_number_, blocStates_, operators);
            transform_operators(transformation_, operators);
            
            jOperators = jsx::array(flavor_number_);
            
            for(unsigned int f = 0; f < operators.size(); ++f) {
                jsx::array jOperator(blocStates_.size());
                for(unsigned int sector = 0; sector < blocStates_.size(); ++sector) {
                    jsx::object jBloc;
                    
                    if(operators[f][sector].target() != -1) {
                        jBloc["target"] = static_cast<std::int64_t>(operators[f][sector].target());
                        operators[f][sector].matrix().b64() = b64;
                        jBloc["matrix"] = std::move(operators[f][sector].matrix());
                    } else
                        jBloc["target"] = jsx::null();
                    
                    jOperator[sector] = std::move(jBloc);
                }
                jOperators[f] = std::move(jOperator);
            }
        };
        
        
        
        /*
        template<class Obs>
        void add_observable(Obs const& obs, jsx::value& jObs, int b64 = true) const {
            
            std::vector<io::rmat> observable;
            
            construct_observable(obs, blocStates_, observable);
            
            transform_observable(transformation_, observable);
            
            jsx::array temp;
            for(auto& bloc : observable) {
                bloc.b64() = b64;
                temp.push_back(bloc);
            }
            jObs = std::move(std::move(temp));
        };
        */
    private:
        int const flavor_number_;
        BlocStates blocStates_;
        std::vector<io::rmat> transformation_;
    };
};

#endif






