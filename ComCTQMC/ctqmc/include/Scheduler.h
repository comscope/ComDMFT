#ifndef SCHEDULER_H
#define SCHEDULER_H

#include <ratio>
#include <chrono>
#include <ctime>
#include <tuple>
#include <random>


namespace mc {
    
    enum class Phase { Initialize, Cycle, Sample, Finalize };
    
    struct Scheduler {
        Scheduler() = delete;
        Scheduler(bool thermalised, Phase phase) : thermalised_(thermalised), phase_(phase) {};
        bool thermalised() const { return thermalised_;};
        Phase& phase() { return phase_;}
        virtual bool done() { return true;};
        virtual ~Scheduler() = default;
    protected:
        bool const thermalised_;
        Phase phase_;
    };
    
    struct TimeScheduler : Scheduler {
        TimeScheduler(std::int64_t duration, bool thermalised, Phase phase) :
        Scheduler(thermalised, phase),
        duration_(60.*duration),
        start_(std::chrono::steady_clock::now()) {
        };
        ~TimeScheduler() = default;
        
        bool done() {
            return !(std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_).count() < duration_);
        };
    private:
        double const duration_;
        std::chrono::steady_clock::time_point const start_;
    };
    
    struct StepsScheduler : Scheduler {
        StepsScheduler(std::int64_t stop, bool thermalised, Phase phase) :
        Scheduler(thermalised, phase),
        stop_(stop),
        steps_(0) {
        };
        ~StepsScheduler() = default;
        
        bool done() {
            return ++steps_ >= stop_;
        };
    private:
        std::int64_t const stop_;
        std::int64_t steps_;
    };
    
}


#endif











