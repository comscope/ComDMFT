#ifndef WEIGHT
#define WEIGHT

#include <vector>

#include "Utilities.h"
#include "Trace.h"
#include "Bath.h"
#include "Dynamic.h"
#include "Config.h"
#include "Data.h"

// Gopfertami config .....

namespace we {
    
    template<class HybFunc>
    struct Weight {
        
        static Weight* Instance(jsx::value const& jParams, da::Data<HybFunc> const& data, int const id) {
            return new Weight(jParams, data, id);
        }
        
        static void Destroy(Weight* instance) {
            delete instance;
        }
   
    private:
        
        Weight() = delete;
        Weight(jsx::value const& jParams, da::Data<HybFunc> const& data, int const id) :
        id_(id),
        trace_(jParams, data.eig(), data.ops()),
        baths_(data.hyb().blocks()),
        dyn_(data.dyn() != nullptr ? new dy::Trace(*data.dyn()) : new dy::Abstract()) {
            //----------------------------------------------------------- Read in configurations ------------------------------------------------------------------
            std::ifstream file(("config_" + std::to_string(id) + ".json").c_str());
            
            if(file) {
                std::cout << "Reading from config file ... " << std::flush;
                
                jsx::value jConfig; jsx::read(file, jConfig);
                
                if(jConfig.size() != baths_.size())
                    throw std::runtime_error("MarkovChain: missmatch in config and bath size.");
                
                for(std::size_t b = 0; b < baths_.size(); ++b) {
                    baths_[b] = std::move(ba::Bath(jConfig(b).object(), data.hyb()));   //emplace_back waer au en moglichkeit ...
                    
                    auto const& opsL = baths_[b].opsL();
                    auto const& opsR = baths_[b].opsR();
                    
                    for(std::size_t i = 0; i < opsL.size(); ++i) {
                        if(!trace_.insert(opsR[i].key(), opsR[i].flavor()))
                            throw std::runtime_error("MarkovChain: key in config file appears twice.");
                        if(!trace_.insert(opsL[i].key(), opsL[i].flavor()))
                            throw std::runtime_error("MarkovChain: key in config file appears twice.");

                        dyn_->insert(opsR[i].key(), opsR[i].flavor()%2);
                        dyn_->insert(opsL[i].key(), opsL[i].flavor()%2);
                    }
                }
                
                std::cout << "Ok" << std::endl;
            }
            
            dyn_->ratio();
            dyn_->accept();
            
            file.close();
        };
        Weight(Weight const&) = delete;
        Weight(Weight&&) = delete;
        Weight& operator=(Weight const&) = delete;
        Weight& operator=(Weight&&) = delete;
        
        ~Weight() {
            jsx::array jConfig(baths_.size());
            for(std::size_t i = 0; i < baths_.size(); ++i) {
                jsx::object jEntry;
                baths_[i].save(jEntry);
                jConfig[i] = std::move(jEntry);
            };
            
            std::ofstream file(("config_" + std::to_string(id_) + ".json").c_str());
            jsx::write(jConfig, file);
            file.close();
        };
        
    public:
        
        tr::Flag prepare() {
            if(tr::Flag::Pending != trace_.surviving())
                throw std::runtime_error("Ma: initial trace is zero");
            
            trace_.prepare(0.);
            
            return tr::Flag::Pending;
        };
        
        tr::Flag decide() {
            if(tr::Flag::Pending == trace_.decide()) return tr::Flag::Pending;
            
            trace_.accept();
            
            return tr::Flag::Accept;
        };
        
        tr::Trace& trace() { return trace_;};
        tr::Trace const& trace() const { return trace_;};
        
        std::vector<ba::Bath>& baths() { return baths_;};
        std::vector<ba::Bath> const& baths() const { return baths_;};
        
        dy::Abstract& dyn() { return *dyn_;};
        dy::Abstract const& dyn() const { return *dyn_;};
        
        int sign() const {
            int sign = trace_.sign();
            for(auto const& bath : baths_)
                sign *= bath.det() > .0 ? 1 : -1;
            return sign;
        };
        
        void clean(da::Data<HybFunc> const& data) {
            dyn_->clean();
            for(auto& bath : baths_)
                bath.clean(data.hyb());
        }
    private:
        int const id_;
        
        tr::Trace trace_;
        std::vector<ba::Bath> baths_;
        std::unique_ptr<dy::Abstract> dyn_;
    };
    
}

#endif
