#ifndef MONTECARLO
#define MONTECARLO

#include <ratio>
#include <chrono>
#include <ctime>
#include <random>

#include "Utilities.h"

namespace mc {
    
    enum class State { Store, UpdateA, UpdateB, Critical, Exit, Init };
    
    struct Control {  // usgseh wie arsch von chue das da ....
        Control(jsx::value const& jParams, int id) :
        id_(id),
        measurement_duration_(60.*jParams("measurement time").int64()),
        state_(State::Init),
        duration_(60.*jParams("thermalisation time").int64()),
        start_(std::chrono::steady_clock::now()),
        thermalized_(false),
        thermalization_steps_(0),
        measurement_steps_(0),
        steps_(&thermalization_steps_) {
            std::cout << std::endl;
            
            std::time_t t = std::time(nullptr);
            std::cout << "Start thermalisation of Markov Chain " << id_ << " at " << std::asctime(std::localtime(&t)) << std::endl;
            
            thermalized();
        };
        
        State& state() { return state_;};
        
        bool thermalized() {
            if(thermalized_) return true;
            if(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() < duration_) return false;
            duration_ = measurement_duration_; start_ = std::chrono::steady_clock::now();
            thermalized_ = true; steps_ = &measurement_steps_;
            
            std::time_t t = std::time(nullptr);
            std::cout << "Start measurements of Markov Chain " << id_ << + " at " << std::asctime(std::localtime(&t)) << std::endl;
            
            return true;
        };
        bool done() { return thermalized_ && !(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() < duration_);};
        std::int64_t& steps() { return *steps_;};
        std::int64_t thermalization_steps() const { return thermalization_steps_;};
        std::int64_t measurement_steps() const { return measurement_steps_;};
    private:
        int const id_;
        double const measurement_duration_;
        
        State state_;
            
        double duration_;
        std::chrono::steady_clock::time_point start_;
            
        bool thermalized_;
        std::int64_t thermalization_steps_, measurement_steps_, *steps_;
        
    };

	
    template<class Comm, class MarkovChain>
    void MonteCarlo(std::vector<MarkovChain*>& markovChains, jsx::value& jSimulation)
    {
        
        ut::RandomNumberGenerator<std::mt19937, std::uniform_real_distribution<double> > urng(std::mt19937(12345), std::uniform_real_distribution<double>(.0, 1.));
        
        jsx::value const& jParams = jSimulation("Parameters");
        
        std::int64_t const sweep = jParams.is("sweep") ? jParams("sweep").int64() : 50;
        std::int64_t const store = sweep*(jParams.is("store") ? jParams("store").int64() : 100);
        std::int64_t const clean = sweep*(jParams.is("clean") ? jParams("clean").int64() : 100);
        
        std::int64_t thermalization_steps = 0;
        std::int64_t measurement_steps = 0;
        
        std::size_t const NStreams = markovChains.size();
        std::vector<int> all_streams(NStreams); std::iota(all_streams.begin(), all_streams.end(), 0);
        
        std::vector<Control*> controls(NStreams, 0);
        for(auto stream : all_streams) controls[stream] = new Control(jParams, markovChains[stream]->id());
        
        std::vector<int> active_streams = all_streams; std::size_t index = 0;
        
        int critical_stream = active_streams[0];
        int critical_steps = -1; //es bitz pingelig
        
        while(active_streams.size()) {
            int const stream = active_streams[index];
            
            if(Comm::is_ready(stream)) {
                MarkovChain& markovChain = *markovChains[stream]; Control& control = *controls[stream];
                
                try {
                    switch (control.state()) {
                        case State::Store:
                            markovChain.store(jSimulation["Measurements"]);
                            control.state() = State::UpdateA;
                            
                        case State::UpdateA:
                            if(!markovChain.udpate()) break;
                            if(++control.steps() % clean == 0) markovChain.clean();
                            if(!control.thermalized()) break;
                            if(++critical_steps >= sweep && stream == critical_stream) {
                                control.state() = State::Critical;
                                critical_steps = 0; break;
                            }
                            
                        case State::UpdateB:
                            if(control.done()) {
                                control.state() = State::Exit;
                                markovChain.sample(); break;
                            }
                            if(control.steps() % sweep == 0) markovChain.sample();
                            control.state() = control.steps() % store == 0 ? State::Store : State::UpdateA;
                            break;
                            
                        case State::Critical:
                            if(!markovChain.csample()) break;
                            critical_stream = urng()*active_streams.size();
                            control.state() = State::UpdateB;
                            break;
                            
                        case State::Exit:
                            active_streams.erase(active_streams.begin() + index);
                            thermalization_steps += control.thermalization_steps();
                            measurement_steps += control.measurement_steps();
                            delete controls[stream]; controls[stream] = 0;
                            break;
                            
                        case State::Init:
                            if(!markovChain.udpate()) break;
                            control.state() = State::UpdateA;
                            break;
                            
                    }
                    Comm::launch();
                }
                
                catch(ut::out_of_memory error) {
                    std::cout << "Markov Chain " << markovChain.id() << " killed." << std::endl;
                    
                    if(stream == critical_stream && active_streams.size()) critical_stream = active_streams[0];
                    active_streams.erase(active_streams.begin() + index); index = 0;
                    all_streams.erase(std::remove(all_streams.begin(), all_streams.end(), stream), all_streams.end());
                    
                    thermalization_steps += control.thermalization_steps();
                    measurement_steps += control.measurement_steps();
                    
                    delete markovChains[stream]; markovChains[stream] = 0;
                    delete controls[stream]; controls[stream] = 0;
                }
            }
        
            if(!(++index < active_streams.size())) index = 0;
        }

        if(!all_streams.size()) throw std::runtime_error("MC: no Markov Chain alive after measurements !");

        for(auto stream : all_streams) markovChains[stream]->store(jSimulation["Measurements"]);

        markovChains[all_streams.back()]->cstore(jSimulation["Measurements"]);
	}
}

#endif











