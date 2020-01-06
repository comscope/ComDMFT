#ifndef MONTECARLO_H
#define MONTECARLO_H

#include <ratio>
#include <chrono>
#include <ctime>
#include <tuple>
#include <random>

#include "Scheduler.h"
#include "MarkovChain.h"
#include "Observables.h"

namespace mc {
    
    struct Simulations {
        template<typename Alloc, typename HybVal>
        Simulations(jsx::value& jParams, data::Data const& data, int mcId, ut::Options<Alloc, HybVal> options) :
        batcher_(new imp::Batcher<Alloc>(8192)),
        init_(new state::Init<Alloc>()),
        state_(jParams, data, mcId, options),
        markovChain_(jParams, mcId) {
        }
        
        imp::itf::Batcher* batcher() { return batcher_.get();};
        state::itf::Init& init() { return *init_;};
        state::State& state() { return state_;};
        MarkovChain& markovChain() { return markovChain_;};
        
        std::unique_ptr<Scheduler>& scheduler() { return scheduler_;};
        std::vector<obs::Observables*>& obs() { return obs_;};
        
    private:
        std::unique_ptr<imp::itf::Batcher> batcher_;
        std::unique_ptr<state::itf::Init> init_;
        state::State state_;
        MarkovChain markovChain_;

        std::unique_ptr<Scheduler> scheduler_;
        std::vector<obs::Observables*> obs_;
    };

    
    template<typename Alloc, typename HybVal>
    void MonteCarlo(jsx::value& jParams, jsx::value& jSimulation)
    {
        data::Data data(jParams, ut::Options<Alloc, HybVal>());

        std::vector<std::unique_ptr<Simulations>> simulations(jSimulation.size());
        for(int stream = 0; stream < jSimulation.size(); ++stream)
            simulations[stream].reset(new Simulations(jParams, data, jSimulation(stream).int64(), ut::Options<Alloc, HybVal>()));

        obs::Observables observablesA(jParams, jParams.is("storeA") ? jParams("storeA").int64() : 100);
        obs::Observables observablesB(jParams, jParams.is("storeB") ? jParams("storeB").int64() :  20);
        
        
        for(int stream = 0; stream < jSimulation.size(); ++stream) {
            simulations[stream]->scheduler().reset(new Scheduler(false, Phase::Initialize));
            addUpdates<Alloc, HybVal>(jParams, data, simulations[stream]->markovChain());
        }
        
        obs::addObservablesA<Alloc, HybVal>(jParams, data, observablesA);
        obs::addObservablesB<Alloc, HybVal>(jParams, data, observablesB);
        
        
        std::size_t const sweepA = jParams.is("sweepA") ? jParams("sweepA").int64() : 50;
        std::size_t const sweepB = jParams.is("sweepB") ? jParams("sweepB").int64() : 250;
        
        std::int64_t thermalisationSteps = 0, measurementSteps = 0, stepsA = 0, stepsB = 0, stream = 0;
        
        jsx::array_t jConfigs; jSimulation = jsx::empty_t();

        while(simulations.size()) {
            auto* batcher = simulations[stream]->batcher();

            if(batcher->is_ready()) {
                auto& sim = *simulations[stream];

                try {
                    switch (sim.scheduler()->phase()) {
                        case Phase::Cycle:
                            if(!sim.markovChain().cycle(data, sim.state(), *batcher)) break;

                            if(sim.scheduler()->done()) {
                                if(sim.scheduler()->thermalised()) {
                                    sim.scheduler()->phase() = Phase::Finalize;
                                } else {
                                    thermalisationSteps += sim.markovChain().steps();
                                    if(jParams.is("measurement steps"))
                                        sim.scheduler().reset(new StepsScheduler(jParams("measurement steps").int64(), true, Phase::Cycle));
                                    else
                                        sim.scheduler().reset(new TimeScheduler(jParams("measurement time").int64(), true, Phase::Cycle));
                                }
                                break;
                            }
                            
                            if(sim.scheduler()->thermalised()) {
                                ++stepsA;
                                if(!observablesA.lock() && stepsA >= sweepA) {
                                    observablesA.lock() = true; sim.obs().push_back(&observablesA);
                                    stepsA = 0; sim.scheduler()->phase() = Phase::Sample;
                                };
                                
                                ++stepsB;
                                if(!observablesB.lock() && stepsB >= sweepB) {
                                    observablesB.lock() = true; sim.obs().push_back(&observablesB);
                                    stepsB = 0; sim.scheduler()->phase() = Phase::Sample;
                                };
                            }

                            break;
                            
                        case Phase::Sample:
                            if(!sim.obs().back()->sample(data, sim.state(), jSimulation["measurements"], *batcher)) break;
                            
                            sim.obs().back()->lock() = false; sim.obs().pop_back();
                            if(!sim.obs().size()) sim.scheduler()->phase() = Phase::Cycle;
                            
                            break;
                        
                        case Phase::Finalize:
                            jConfigs.push_back(sim.state().config().json());

                            measurementSteps += sim.markovChain().steps();
                            simulations.erase(simulations.begin() + stream);
                            batcher = nullptr;

                            break;
                            
                        case Phase::Initialize:
                            if(!sim.init().apply(data, sim.state(), *batcher)) break;

                            if(jParams.is("thermalisation steps"))
                                sim.scheduler().reset(new StepsScheduler(jParams("thermalisation steps").int64(), false, Phase::Cycle));
                            else
                                sim.scheduler().reset(new TimeScheduler(jParams("thermalisation time").int64(), false, Phase::Cycle));

                            break;
                    }

                    if(batcher != nullptr) batcher->launch();
                }
                
                catch(ut::out_of_memory error) {
                    if(sim.scheduler()->phase() == Phase::Sample) throw std::runtime_error("MC: Fatal error, out of memory while sampling");
                    
                    std::cout << "MC: Markov Chain gets killed." << std::endl;

                    measurementSteps += sim.markovChain().steps();
                    simulations.erase(simulations.begin() + stream);
                    
                    if(!simulations.size()) throw std::runtime_error("MC: all Markov Chain killed ! Himmelhergottsakramentzifixhallelujaleckstmiamarscheissglumpvereckts !");
                }
            }
            if(!(++stream < simulations.size())) stream = 0;
        }

        observablesA.finalize(data, jSimulation["measurements"]);
        observablesB.finalize(data, jSimulation["measurements"]);
        jSimulation["configs"] = std::move(jConfigs);
        
        std::int64_t numberOfMarkovChains = jSimulation["configs"].size();
        mpi::reduce<mpi::op::sum>(numberOfMarkovChains, mpi::master);
        mpi::reduce<mpi::op::sum>(thermalisationSteps, mpi::master);
        mpi::reduce<mpi::op::sum>(measurementSteps, mpi::master);

        jSimulation["info"]["number of mpi processes"] = jsx::int64_t(mpi::number_of_workers());
        jSimulation["info"]["number of markov chains"] = numberOfMarkovChains;
        jSimulation["info"]["thermalization steps"] = thermalisationSteps;
        jSimulation["info"]["measurement steps"] = measurementSteps;
    }
}


#endif











