#ifndef CTQMC_INCLUDE_MONTECARLO_H
#define CTQMC_INCLUDE_MONTECARLO_H

#include <ratio>
#include <chrono>
#include <ctime>
#include <tuple>
#include <random>

#include "Params.h"

#include "markovchain/Scheduler.h"
#include "markovchain/MarkovChain.h"

#include "updates/Setup.h"

#include "observables/Observables.h"
#include "observables/Setup.h"

#include "../../include/io/Tag.h"
#include "../../include/measurements/Error.h"
#include "../../evalsim/Evalsim.h"

namespace mc {


    template<typename Mode, typename Value>
    void montecarlo(jsx::value jParams, jsx::value& jSimulation)
    {
        params::complete_impurity<Value>(jParams);
        
        data::Data<Value> data(jParams, Mode());
        data::setup_data<Mode>(jParams, data);
        
        obs::Observables<Value> observables;
        obs::setup_obs<Mode>(jParams, data, observables);
        
        mch::WangLandau<Value> wangLandau(jParams, data);
        
        std::vector<std::tuple<
        std::unique_ptr<imp::itf::Batcher<Value>>, //0
        std::unique_ptr<state::State<Value>     >, //1
        std::unique_ptr<mch::MarkovChain<Value> >, //2
        std::unique_ptr<mch::Scheduler          >  //3
        >> simulations;
        
        for(int stream = 0; stream < jSimulation.size(); ++stream) {
            
            simulations.emplace_back(
            std::unique_ptr<imp::itf::Batcher<Value>>(new imp::Batcher<Mode, Value>(8192)),
            std::unique_ptr<state::State<Value>     >(new state::State<Value>(jParams, data, jSimulation(stream)("config"), Mode())),
            std::unique_ptr<mch::MarkovChain<Value> >(new mch::MarkovChain<Value>(jParams, jSimulation(stream)("id").int64(), Mode())),
            std::unique_ptr<mch::Scheduler          >(new mch::Scheduler(false, mch::Phase::Initialize))
            );

            upd::setup_updates<Mode>(jParams, data, *std::get<1>(simulations.back()), *std::get<2>(simulations.back()));
        }

        jSimulation["configs"] = jsx::array_t();
        
        
        std::int64_t thermSteps = 0, measSteps = 0, stream = 0;
        if(jParams.is("restart") and jParams("restart").boolean()){
            meas::restart(jParams, jParams("measurements"), jSimulation["measurements"]); jParams["measurements"] = jsx::empty_t();
        }
            //I kind of want to free the memory from jParams[measurements] earlier. However, code complains because of jSimulation
        
        while(simulations.size()) {
            auto* batcher = std::get<0>(simulations[stream]).get();

            if(batcher->is_ready()) {
                auto& state       = std::get<1>(simulations[stream]);
                auto& markovChain = std::get<2>(simulations[stream]);
                auto& scheduler   = std::get<3>(simulations[stream]);

                try {
                    switch (scheduler->phase()) {
                        case mch::Phase::Step:
                            if(!markovChain->cycle(wangLandau, data, *state, *batcher)) break;
                                         
                            if(scheduler->done()) {
                                if(scheduler->thermalised()) {
                                    scheduler->phase() = mch::Phase::Finalize;
                                } else {
                                    if(jParams.is("measurement steps"))
                                        scheduler.reset(new mch::StepsScheduler(jParams("measurement steps").int64(), true, mch::Phase::Step));
                                    else
                                        scheduler.reset(new mch::TimeScheduler(jParams("measurement time").int64(), true, mch::Phase::Step));
                                    wangLandau.thermalised();
                                }
                                break;
                            }
                            
                            if(scheduler->thermalised()) {
                                ++measSteps;
                                
                                if(observables[state->worm().index()]->sample(data, *state))
                                    scheduler->phase() = mch::Phase::Sample;
                            } else
                                ++thermSteps;

                            break;
                            
                        case mch::Phase::Sample:
                            if(!observables[state->worm().index()]->cycle(data, *state, jSimulation["measurements"], *batcher)) break;
                            
                            scheduler->phase() = mch::Phase::Step;
                            
                            break;
                        
                        case mch::Phase::Finalize:
                            jSimulation["configs"].array().push_back(state->json());

                            simulations.erase(simulations.begin() + stream);
                            batcher = nullptr;

                            break;
                            
                        case mch::Phase::Initialize:
                            if(!markovChain->init(data, *state, *batcher)) break;

                            if(jParams.is("thermalisation steps"))
                                scheduler.reset(new mch::StepsScheduler(jParams("thermalisation steps").int64(), false, mch::Phase::Step));
                            else
                                scheduler.reset(new mch::TimeScheduler(jParams("thermalisation time").int64(), false, mch::Phase::Step));
                            
                            break;
                    }

                    if(batcher != nullptr) batcher->launch();
                }
                
                catch(ut::out_of_memory error) {
                    if(std::get<3>(simulations[stream])->phase() == mch::Phase::Sample)
                        throw std::runtime_error("MC: Fatal error, out of memory while sampling");
                    
                    std::cout << "MC: Markov Chain gets killed." << std::endl;

                    simulations.erase(simulations.begin() + stream);
                    
                    if(!simulations.size())
                        throw std::runtime_error("MC: all Markov Chain killed !");
                }
            }
            if(!(++stream < simulations.size())) stream = 0;
        }
        
        for(std::size_t space = 0; space < cfg::Worm::size(); ++space)
            if(observables[space] != nullptr)
                observables[space]->finalize(data, jSimulation["measurements"]);
        
        jSimulation["Wang Landau"] = wangLandau.json();
        
        std::int64_t numberOfMarkovChains = jSimulation["configs"].size();
        
        mpi::reduce<mpi::op::sum>(numberOfMarkovChains, mpi::master);
        mpi::reduce<mpi::op::sum>(thermSteps,           mpi::master);
        mpi::reduce<mpi::op::sum>(measSteps,            mpi::master);

        jSimulation["info"] = jsx::object_t{
            { "number of mpi processes", mpi::number_of_workers() },
            { "number of markov chains", numberOfMarkovChains },
            { "thermalization steps",    thermSteps },
            { "measurement steps",       measSteps }
        };

    }
    
    
    
    
    template<typename Value>
    void statistics(jsx::value jParams, jsx::value& jSimulation) {
        if(mpi::number_of_workers() > 1 && jParams.is("error") && jParams("error").string() != "none") {
            jsx::value jMeasurements = std::move(jSimulation("measurements"));
            
            meas::reduce(jSimulation("measurements"), jMeasurements, jSimulation("Wang Landau"), meas::All(), true);

            if(jParams("error").string() == "parallel") {
                
                meas::reduce(jMeasurements, jMeasurements, jSimulation("Wang Landau"), meas::Jackknife(), false);
                jSimulation["error"] = evalsim::evalsim<Value>(jParams, jMeasurements);
                meas::error(jSimulation("error"), meas::Jackknife());
                
            } else if(jParams("error").string() == "serial") {
                
                meas::reduce(jMeasurements, jMeasurements, jSimulation("Wang Landau"), meas::Jackknife(), true);
                jSimulation["resample"] = std::move(jMeasurements);
                io::to_tagged_json(jSimulation("resample"));
                
            } else
                throw std::runtime_error("mc::statistics: invalid error option " + jParams("error").string());
            
        } else
            meas::reduce(jSimulation("measurements"), jSimulation("measurements"), jSimulation("Wang Landau"), meas::All(), true);
        
        io::to_tagged_json(jSimulation("measurements"));
    }
}


#endif











