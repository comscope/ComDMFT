#ifndef INCLUDE_ATOMIC_GENERATE_H
#define INCLUDE_ATOMIC_GENERATE_H

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

#include "Tensor.h"
#include "../JsonX.h"
#include "../io/Vector.h"
#include "../io/Matrix.h"
#include "../linalg/LinAlg.h"
#include "../mpi/Utilities.h"

namespace ga {
    //kack isch veralted ....
    
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
    typedef std::vector<States> BlockStates;
    
    struct Join {
        Join(std::size_t size) : labels_(size) {
            std::iota(labels_.begin(), labels_.end(), 0);
        };
        void join(State state1, State state2) {
            state1 = find_representative(state1);
            state2 = find_representative(state2);
            state1 < state2 ? labels_[state2] = state1 : labels_[state1] = state2;
        };
        std::size_t clean() {
            for(auto & label : labels_) label = find_representative(label);
            std::map<State, State> new_label;
            for(auto & label : labels_) {
                if(!new_label.count(label)) {
                   State temp = new_label.size();
                   new_label[label] = temp;
                }
                label = new_label[label];
            }
            return new_label.size();
        };
        State label(State state) {
            return labels_[state];
        };
    private:
        std::vector<State> labels_;
        
        State find_representative(State label) {
            while(label != labels_[label]) label = labels_[label];
            return label;
        };
    };
    
    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    
    template<typename Value>
    void find_blocks(Tensor<Value> const& hloc, BlockStates& blockStates)
    {
        Join pre_block_labels(1 << hloc.N());
      
        for(State state = 0; state < (1 << hloc.N()); ++state) {
            FlavorState const stateJ(state);
            
            for(int fDagg = 0; fDagg < hloc.N(); ++fDagg)
                for(int f = 0; f < hloc.N(); ++f)
                    if(std::abs(hloc.t(fDagg, f)) > 1.e-14) {
                        FlavorState const stateI = psiDagg(fDagg, psi(f, stateJ));
                        
                        if(stateI.sign() != 0)
                            pre_block_labels.join(stateI.state(), stateJ.state());
                    }
            
            for(int f1Dagg = 0; f1Dagg < hloc.N(); ++f1Dagg)
                for(int f2Dagg = 0; f2Dagg < hloc.N(); ++f2Dagg)
                    for(int f1 = 0; f1 < hloc.N(); ++f1)
                        for(int f2 = 0; f2 < hloc.N(); ++f2)
                            if(std::abs(hloc.V(f1Dagg, f2Dagg, f1, f2)) > 1.e-14) {
                                FlavorState const stateI = psiDagg(f1Dagg, psiDagg(f2Dagg, psi(f1, psi(f2, stateJ))));
                                
                                if(stateI.sign() != 0)
                                    pre_block_labels.join(stateI.state(), stateJ.state());
                            }
        }

        Join block_labels(pre_block_labels.clean());
        
        for(int f = 0; f < hloc.N(); ++f)
            for(int type = 0; type < 2; ++type) {
                std::map<State, State> mapping;
                for(State state = 0; state < (1 << hloc.N()); ++state) {
                    FlavorState const stateJ(state);
                    FlavorState const stateI = type ? psiDagg(f, stateJ) : psi(f, stateJ);
                    if(stateI.sign() != 0) {
                        int const blockJ = pre_block_labels.label(stateJ.state());
                        int const blockI = pre_block_labels.label(stateI.state());
                        
                        if(mapping.count(blockJ))
                            block_labels.join(blockI, mapping[blockJ]);
                        else
                            mapping[blockJ] = blockI;
                    }
                }
            }
       
                                                  
       blockStates.resize(block_labels.clean());
       for(State state = 0; state < (1 << hloc.N()); ++state) 
            blockStates[block_labels.label(pre_block_labels.label(state))].push_back(state);
    };

    
    enum class Order { alternating, normal };
    
    
    template<Order order, typename Value>
    jsx::value get_observable(Tensor<Value> const& tensor,
                              BlockStates const& blockStates,
                              bool const throw_error = true) //if false, instead return jsx::empty on error
    {
        jsx::value jObservable = jsx::array_t(blockStates.size());

        State block_label = 0;
        for(auto const& blockState : blockStates) {
            jObservable[block_label]["target"] = jsx::int64_t(block_label);
            io::Matrix<Value> matrix(blockState.size(), blockState.size());
            
            State state_indexJ = 0;
            for(auto const & state : blockState) {
                FlavorState const stateJ(state);
                
                for(int f1 = 0; f1 < tensor.N(); ++f1)
                    for(int f2 = 0; f2 < tensor.N(); ++f2)
                        if(std::abs(tensor.t(f1, f2)) > 1.e-14) {
                            FlavorState const stateI = psiDagg(f1, psi(f2, stateJ));
                            
                            if(stateI.sign() != 0) {
                                auto it = std::find(blockState.begin(), blockState.end(), stateI.state());
                                
                                if(it == blockState.end()){
                                    if (throw_error) throw std::runtime_error("Something is wrong with the partitioning of the states");
                                    else {mpi::cout << " not a good observable; removing from list ... " ; return jsx::empty();}
                                }
                                
                                
                                matrix(it - blockState.begin(), state_indexJ) += tensor.t(f1, f2)*static_cast<double>(stateI.sign());
                            }
                        }
                
                for(int f1 = 0; f1 < tensor.N(); ++f1)
                    for(int f2 = 0; f2 < tensor.N(); ++f2)
                        for(int f3 = 0; f3 < tensor.N(); ++f3)
                            for(int f4 = 0; f4 < tensor.N(); ++f4)
                                if(std::abs(tensor.V(f1, f2, f3, f4)) > 1.e-14) {
                                    FlavorState const stateI = order == Order::normal ? psiDagg(f1, psiDagg(f2, psi(f3, psi(f4, stateJ)))) : psiDagg(f1, psi(f2, psiDagg(f3, psi(f4, stateJ))));
                                    
                                    if(stateI.sign() != 0) {
                                        auto it = std::find(blockState.begin(), blockState.end(), stateI.state());
                                        
                                        if(it == blockState.end()){
                                            if (throw_error) throw std::runtime_error("Something is wrong with the partitioning of the states");
                                            else {mpi::cout << " not a good observable; removing from list ... " ; return jsx::empty();}
                                        }
                                        
                                        matrix(it - blockState.begin(), state_indexJ) += tensor.V(f1, f2, f3, f4)*static_cast<double>(stateI.sign());
                                    }
                                }
                ++state_indexJ;
            }
            jObservable[block_label]["matrix"] = std::move(matrix);
        
            ++block_label;
        }
        
        return jObservable;
    };

    
    template<typename Value>
    jsx::value diagonalise(jsx::value& jHamiltonian)
    {
        jsx::value jEigenValues = jsx::array_t(jHamiltonian.size());

        for(unsigned int sector = 0; sector < jHamiltonian.size(); ++sector) {
            jEigenValues(sector) = io::rvec(jsx::at<io::Matrix<Value>>(jHamiltonian(sector)("matrix")).I());
            linalg::eig('V', 'U', jsx::at<io::Matrix<Value>>(jHamiltonian(sector)("matrix")), jsx::at<io::rvec>(jEigenValues[sector]));
        }
        
        return jEigenValues;
    };
    
    
    template<typename Value>
    void transform(jsx::value const& jTransformation,
                   jsx::value& jOperator)
    {
        for(std::size_t start = 0; start < jTransformation.size(); ++start)
            if(!jOperator(start)("target").is<jsx::null_t>()) {
                auto const target = jOperator(start)("target").int64();
                io::Matrix<Value> buffer(jsx::at<io::Matrix<Value>>(jOperator(start)("matrix")).I(), jsx::at<io::Matrix<Value>>(jOperator(start)("matrix")).J());
                
                linalg::mult<Value>('n', 'n', 1., jsx::at<io::Matrix<Value>>(jOperator(start)("matrix")), jsx::at<io::Matrix<Value>>(jTransformation(start)("matrix")), .0, buffer);
                linalg::mult<Value>('c', 'n', 1., jsx::at<io::Matrix<Value>>(jTransformation(target)("matrix")), buffer, .0, jsx::at<io::Matrix<Value>>(jOperator(start)("matrix")));
            }
    };
    
    
    template<typename Value>
    jsx::value get_annihilation_operators(int N,
                                          BlockStates const& blockStates)
    {
        jsx::value jOperators = jsx::array_t(N, jsx::array_t(blockStates.size(), jsx::object_t{{"target", jsx::null_t()}}));

        std::vector<State> block_labels(1 << N); State block_label = 0;
        for(auto const& blockState : blockStates) {
            for(auto const& state : blockState) block_labels[state] = block_label;
            ++block_label;
        }
        
        State block_labelJ = 0;
        for(auto const& blockStateJ : blockStates) {
            State state_indexJ = 0;
            for(auto const& state : blockStateJ) {
                FlavorState const stateJ(state);
                for(int f = 0; f < N; ++f) {
                    FlavorState const stateI = psi(f, stateJ);
                    if(stateI.sign() != 0) {
                        State const block_labelI = block_labels[stateI.state()];
                        auto const& blockStateI = blockStates[block_labelI];
                        State const state_indexI = std::find(blockStateI.begin(), blockStateI.end(), stateI.state()) - blockStateI.begin();
                        
                        if(jOperators(f)(block_labelJ)("target").is<jsx::null_t>()) {
                            jOperators(f)(block_labelJ)("target") = jsx::int64_t(block_labelI);
                            jOperators(f)(block_labelJ)["matrix"] = io::Matrix<Value>(blockStateI.size(), blockStateJ.size());
                        } else if(jOperators(f)(block_labelJ)("target").int64() != static_cast<jsx::int64_t>(block_labelI))
                            throw std::runtime_error("Fatal error: Operator maps one sector to two or more sectors.");
                        
                        jsx::at<io::Matrix<Value>>(jOperators(f)(block_labelJ)("matrix"))(state_indexI, state_indexJ) += stateI.sign();
                    }
                }
                ++state_indexJ;
            }
            ++block_labelJ;
        }
        
        return jOperators;
    };

    io::rvec get_sector_qn(BlockStates const& blockStates,
                             std::vector<double> const& qn,
                           bool const throw_error = true)
    {
        io::rvec Qn;
        double const eps = 1E-4;
        
        for(auto const& states : blockStates) {
            double value = .0;
            
            for(std::size_t f = 0; f < qn.size(); ++f)
                if(FlavorState(states.front()).test(f)) value += qn[f];
            
            for(auto const state : states) {
                double temp = .0;
                
                for(std::size_t f = 0; f < qn.size(); ++f)
                    if(FlavorState(state).test(f)) temp += qn[f];
                
                if(std::abs(value  - temp) > eps){
                    if (throw_error) throw std::runtime_error("Problem with quantum numbers.");
                    else {mpi::cout << " not a good qn; removing from list " ; return io::rvec();}
                }
                    
            }
            
            Qn.push_back(value);
        }
        
        return Qn;
    };
    
    //-----------------------------------------------------------------------------------------------------
    //-----------------------------------------------------------------------------------------------------
    
    template<typename Value>
    jsx::value construct_hloc(jsx::value jTensors, bool b64 = true)
    {
        jsx::value jHloc; Tensor<Value> hloc(jTensors);
        
        BlockStates blockStates;
        mpi::cout << "Number of invariant subspaces: " << std::flush;
        find_blocks(hloc, blockStates);
        mpi::cout << blockStates.size() << std::endl;
        
        std::size_t maxDim = 0;
        for(auto const& states : blockStates) maxDim = std::max(states.size(), maxDim);
        mpi::cout << "Dimension of the biggest subspace: " << maxDim << std::endl;

        jHloc["transformation"] = get_observable<Order::normal>(hloc, blockStates);
        jHloc["eigen values"] = diagonalise<Value>(jHloc("transformation"));

        jHloc["interaction"] = get_observable<Order::normal>(Tensor<Value>(hloc, typename Tensor<Value>::Interaction()), blockStates);
        transform<Value>(jHloc("transformation"), jHloc("interaction"));
        
        jHloc["filling"] = get_sector_qn(blockStates, std::vector<double>(hloc.N(), 1.));
        
        io::Matrix<Value> one_body(hloc.N(), hloc.N()); one_body.b64() = b64;
        for(int fDagger = 0; fDagger < hloc.N(); ++fDagger)
            for(int f = 0; f < hloc.N(); ++f)
                one_body(fDagger, f) = hloc.t(fDagger, f);
        jHloc["one body"] = std::move(one_body);
        jHloc["two body"] = jsx::at<io::Vector<Value>>(jTensors("two body"));

        jsx::array_t jBlockStates;
        for(auto const& states : blockStates) {
            jsx::array_t jStates;
            for(auto const& state : states)
                jStates.push_back(jsx::int64_t(state));
            jBlockStates.push_back(std::move(jStates));
        }
        jHloc["block states"] = std::move(jBlockStates);
        
        jsx::at<io::rvec>(jHloc("filling")).b64() = b64;
        for(auto& jBlock : jHloc("eigen values").array())   jsx::at<io::rvec>(jBlock).b64() = b64;
        for(auto& jBlock : jHloc("transformation").array()) jsx::at<io::Matrix<Value>>(jBlock("matrix")).b64() = b64;
        for(auto& jBlock : jHloc("interaction").array())    jsx::at<io::Matrix<Value>>(jBlock("matrix")).b64() = b64;
        
        return jHloc;
    };
    
    
    //-----------------------------------------------------------------------------------------------------
    
    template<typename Value>
    jsx::value read_hloc(std::string const name)
    {
        jsx::value jHloc = mpi::read(name);
        
        if(jsx::at<io::Matrix<Value>>(jHloc("one body")).I() != jsx::at<io::Matrix<Value>>(jHloc("one body")).J())
            throw std::runtime_error("ga::real_hloc: one body matrix not square");
        
        int const sectorNumber = jHloc("block states").size();
        
        if(jHloc("eigen values").size() != sectorNumber)
            throw std::runtime_error("ga::sanity_check: eigenvalues have wrong sector number");
        
        if(jHloc("transformation").size() != sectorNumber)
            throw std::runtime_error("ga::sanity_check: transformation has wrong sector number");
        
        if(jsx::at<io::rvec>(jHloc("filling")).size() != sectorNumber)
            throw std::runtime_error("ga::sanity_check: transformation has wrong sector number");

        
        for(std::size_t sector = 0; sector < sectorNumber; ++sector) {
            auto const& eigenvalues = jsx::at<io::rvec>(jHloc("eigen values")(sector));
            
            if(jHloc("block states")(sector).size() != eigenvalues.size())
                throw std::runtime_error("ga::sanity_check: eigenvalues and block states not compatible");
            
            auto const& transformation = jsx::at<io::Matrix<Value>>(jHloc("transformation")(sector)("matrix"));
            
            if(jHloc("block states")(sector).size() != static_cast<std::size_t>(transformation.I()) ||
               jHloc("block states")(sector).size() != static_cast<std::size_t>(transformation.J()))
                throw std::runtime_error("ga::sanity_check: transformation and block states not compatible");
        }
        
        return jHloc;
    };
    
    
    BlockStates get_block_states(jsx::value const& jHloc)
    {
        BlockStates blockStates;
        
        for(auto const& jStates : jHloc("block states").array()) {
            States states;
            for(auto const& jState : jStates.array())
                states.push_back(jState.int64());
            blockStates.push_back(std::move(states));
        }
        
        return blockStates;
    };
    
    
    //-----------------------------------------------------------------------------------------------------
    
    io::rvec construct_sector_qn(jsx::value const& jHloc, jsx::value jqn, bool throw_error = true)
    {
        return get_sector_qn(get_block_states(jHloc), jsx::at<io::rvec>(jqn), throw_error);
    };
    
    template<typename Value>
    jsx::value construct_annihilation_operators(jsx::value const& jHloc)
    {
        jsx::value jOperators = get_annihilation_operators<Value>(jsx::at<io::Matrix<Value>>(jHloc("one body")).I(), get_block_states(jHloc));
        
        for(auto& jOperator : jOperators.array()) transform<Value>(jHloc("transformation"), jOperator);
        
        return jOperators;
    };

    template<typename Value>
    jsx::value construct_observable(jsx::value const& jHloc, jsx::value const& jTensors, bool const throw_error = true)
    {
        BlockStates blockStates = get_block_states(jHloc);
        
        jsx::value jObservable = get_observable<Order::alternating>(Tensor<Value>(jTensors), blockStates, throw_error);
        
        if (!throw_error and jObservable.is<jsx::empty_t>()) return jObservable;
        
        transform<Value>(jHloc("transformation"), jObservable);
        
        return jObservable;
    };
    
    template<typename Value>
    jsx::value construct_occupation_states(jsx::value const& jHloc)
    {
        jsx::array_t jOccupationStates;
        
        BlockStates blockStates = get_block_states(jHloc);
        for(auto const& states : blockStates)
            for(auto const state : states) {
                io::rvec temp;
                for(int f = 0; f < jsx::at<io::Matrix<Value>>(jHloc("one body")).I(); ++f)
                    temp.push_back(FlavorState(state).test(f) ? 1. : .0);
                jOccupationStates.push_back(temp);
            }
        
        return jOccupationStates;
    };

};

#endif







