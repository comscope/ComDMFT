#ifndef CTQMC_INCLUDE_UPDATES_SETUP_H
#define CTQMC_INCLUDE_UPDATES_SETUP_H

#include <vector>
#include <array>

#include "expansion/Setup.h"

#include "worm/InsertOps.h"
#include "worm/RemoveOps.h"
#include "worm/Reconnect.h"

#include "../markovchain/Update.h"
#include "../markovchain/MarkovChain.h"


namespace upd {
    
    
    template<typename Update, typename Mode, typename Value, typename...Args>
    mch::unique_update_ptr<Value> make_update(double prob, Args&&... args) {
        return mch::unique_update_ptr<Value>(new Generic<Update, Mode, Value>(prob, std::forward<Args>(args)...));
    }
    
    
    template<typename Mode, typename Value>
    void setup_updates(jsx::value const& jParams, data::Data<Value>& data, state::State<Value> const& state, mch::MarkovChain<Value>& markovChain) {
        
        using namespace cfg;  using namespace worm;
        
        mpi::cout << "Begin setting updates" << std::endl;

        
        expansion::setup_updates<partition::Worm, Mode>(jParams, data, state, markovChain);
        
        
        if(jParams.is(green::Worm::name())) {
            
            expansion::setup_updates<green::Worm, Mode>(jParams, data, state, markovChain);
            
            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<op, opDagg>, green::Worm, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<green::Worm, partition::Worm, std::tuple<op, opDagg>, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
            
            markovChain.add(make_update<Reconnect<green::Worm, 0>, Mode, Value>(.5));
            markovChain.add(make_update<Reconnect<green::Worm, 1>, Mode, Value>(.5));
            
        }
    
        
        if(jParams.is(green_impr::Worm::name())) {
            
            expansion::setup_updates<green_impr::Worm, Mode>(jParams, data, state, markovChain);
            
            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBulla, opDagg>, green_impr::Worm, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<green_impr::Worm, partition::Worm, std::tuple<opBulla, opDagg>, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
            
            markovChain.add(make_update<Reconnect<green_impr::Worm, 1>, Mode, Value>(1.));
            
        }
        
        
        if(jParams.is(green_imprsum::Worm::name())) {
            
            expansion::setup_updates<green_imprsum::Worm, Mode>(jParams, data, state, markovChain);
            
            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBullaSum, opDagg>, green_imprsum::Worm, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<green_imprsum::Worm, partition::Worm, std::tuple<opBullaSum, opDagg>, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<green_imprsum::Worm, 1>, Mode, Value>(1.));
            
        }
        
        
        if(jParams.is(vertex::Worm::name())) {
            
            expansion::setup_updates<vertex::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<op, opDagg, op, opDagg>, vertex::Worm, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<vertex::Worm, partition::Worm, std::tuple<op, opDagg, op, opDagg>, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<vertex::Worm, 0>, Mode, Value>(.25));
            markovChain.add(make_update<Reconnect<vertex::Worm, 1>, Mode, Value>(.25));
            markovChain.add(make_update<Reconnect<vertex::Worm, 2>, Mode, Value>(.25));
            markovChain.add(make_update<Reconnect<vertex::Worm, 3>, Mode, Value>(.25));
            
            
            if(jParams.is(green::Worm::name())) {

                markovChain.add(make_update<InsertOps<green::Worm, std::tuple<op, opDagg>, vertex::Worm, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<vertex::Worm, green::Worm, std::tuple<op, opDagg>, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data));

                markovChain.add(make_update<InsertOps<green::Worm, std::tuple<op, opDagg>, vertex::Worm, ut::sequence<2, 3, 0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<vertex::Worm, green::Worm, std::tuple<op, opDagg>, ut::sequence<2, 3, 0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }

        }
        
        
        if(jParams.is(vertex_impr::Worm::name())) {
            
            expansion::setup_updates<vertex_impr::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBulla, opDagg, op, opDagg>, vertex_impr::Worm, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<vertex_impr::Worm, partition::Worm, std::tuple<opBulla, opDagg, op, opDagg>, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<vertex_impr::Worm, 1>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<vertex_impr::Worm, 2>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<vertex_impr::Worm, 3>, Mode, Value>(1./3.));
            
            
            if(jParams.is(green::Worm::name())) {

                markovChain.add(make_update<InsertOps<green::Worm, std::tuple<opBulla, opDagg>, vertex_impr::Worm, ut::sequence<2, 3, 0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<vertex_impr::Worm, green::Worm, std::tuple<opBulla, opDagg>, ut::sequence<2, 3, 0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }
            
            
            if(jParams.is(green_impr::Worm::name())) {

                markovChain.add(make_update<InsertOps<green_impr::Worm, std::tuple<op, opDagg>, vertex_impr::Worm, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<vertex_impr::Worm, green_impr::Worm, std::tuple<op, opDagg>, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }
            
        }
        
        
        if(jParams.is(vertex_imprsum::Worm::name())) {
            
            expansion::setup_updates<vertex_imprsum::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBullaSum, opDagg, op, opDagg>, vertex_imprsum::Worm, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<vertex_imprsum::Worm, partition::Worm, std::tuple<opBullaSum, opDagg, op, opDagg>, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<vertex_imprsum::Worm, 1>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<vertex_imprsum::Worm, 2>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<vertex_imprsum::Worm, 3>, Mode, Value>(1./3.));
            
            
            if(jParams.is(green::Worm::name())) {

                markovChain.add(make_update<InsertOps<green::Worm, std::tuple<opBullaSum, opDagg>, vertex_imprsum::Worm, ut::sequence<2, 3, 0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<vertex_imprsum::Worm, green::Worm, std::tuple<opBullaSum, opDagg>, ut::sequence<2, 3, 0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }
            
            
            if(jParams.is(green_imprsum::Worm::name())) {

                markovChain.add(make_update<InsertOps<green_imprsum::Worm, std::tuple<op, opDagg>, vertex_imprsum::Worm, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<vertex_imprsum::Worm, green_imprsum::Worm, std::tuple<op, opDagg>, ut::sequence<0, 1, 2, 3>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }
            
        }
        
        
        if(jParams.is(susc_ph::Worm::name())) {
            
            expansion::setup_updates<susc_ph::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<bilinearPH, bilinearPH>, susc_ph::Worm, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<susc_ph::Worm, partition::Worm, std::tuple<bilinearPH, bilinearPH>, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
            
            markovChain.add(make_update<Reconnect<susc_ph::Worm, 0>, Mode, Value>(0.5));
            markovChain.add(make_update<Reconnect<susc_ph::Worm, 1>, Mode, Value>(0.5));
   
        }
        
        
        if(jParams.is(susc_pp::Worm::name())) {
            
            expansion::setup_updates<susc_pp::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<bilinearHH, bilinearPP>, susc_pp::Worm, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<susc_pp::Worm, partition::Worm, std::tuple<bilinearHH, bilinearPP>, ut::sequence<0, 1>>, Mode, Value>(1., jsx::empty_t(), data));
            
            markovChain.add(make_update<Reconnect<susc_pp::Worm, 0>, Mode, Value>(0.5));
            markovChain.add(make_update<Reconnect<susc_pp::Worm, 1>, Mode, Value>(0.5));
            
        }
        
        
        if(jParams.is(hedin_ph::Worm::name())) {
            
            expansion::setup_updates<hedin_ph::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<op, opDagg, bilinearPH>, hedin_ph::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<hedin_ph::Worm, partition::Worm, std::tuple<op, opDagg, bilinearPH>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<hedin_ph::Worm, 0>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<hedin_ph::Worm, 1>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<hedin_ph::Worm, 2>, Mode, Value>(1./3.));
            
            
            if(jParams.is(green::Worm::name())) {

                markovChain.add(make_update<InsertOps<green::Worm, std::tuple<bilinearPH>, hedin_ph::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<hedin_ph::Worm, green::Worm, std::tuple<bilinearPH>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            }
            
        }
        
        
        if(jParams.is(hedin_ph_impr::Worm::name())) {
            
            expansion::setup_updates<hedin_ph_impr::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBulla, opDagg, bilinearPH>, hedin_ph_impr::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<hedin_ph_impr::Worm, partition::Worm, std::tuple<opBulla, opDagg, bilinearPH>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<hedin_ph_impr::Worm, 1>, Mode, Value>(.5));
            markovChain.add(make_update<Reconnect<hedin_ph_impr::Worm, 2>, Mode, Value>(.5));
            
            if(jParams.is(green_impr::Worm::name())) {

                markovChain.add(make_update<InsertOps<green_impr::Worm, std::tuple<bilinearPH>, hedin_ph_impr::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<hedin_ph_impr::Worm, green_impr::Worm, std::tuple<bilinearPH>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }
            
        }
        
        
        if(jParams.is(hedin_ph_imprsum::Worm::name())) {
            
            expansion::setup_updates<hedin_ph_imprsum::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBullaSum, opDagg, bilinearPH>, hedin_ph_imprsum::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<hedin_ph_imprsum::Worm, partition::Worm, std::tuple<opBullaSum, opDagg, bilinearPH>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<hedin_ph_imprsum::Worm, 1>, Mode, Value>(.5));
            markovChain.add(make_update<Reconnect<hedin_ph_imprsum::Worm, 2>, Mode, Value>(.5));
            
            if(jParams.is(green_imprsum::Worm::name())) {

                markovChain.add(make_update<InsertOps<green_imprsum::Worm, std::tuple<bilinearPH>, hedin_ph_imprsum::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                                make_update<RemoveOps<hedin_ph_imprsum::Worm, green_imprsum::Worm, std::tuple<bilinearPH>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));
                
            }
            
        }
        
        
        if(jParams.is(hedin_pp::Worm::name())) {
            
            expansion::setup_updates<hedin_pp::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<op, op, bilinearPP>, hedin_pp::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<hedin_pp::Worm, partition::Worm, std::tuple<op, op, bilinearPP>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<hedin_pp::Worm, 0>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<hedin_pp::Worm, 1>, Mode, Value>(1./3.));
            markovChain.add(make_update<Reconnect<hedin_pp::Worm, 2>, Mode, Value>(1./3.));
            
        }
        
        
        if(jParams.is(hedin_pp_impr::Worm::name())) {
            
            expansion::setup_updates<hedin_pp_impr::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBulla, op, bilinearPP>, hedin_pp_impr::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<hedin_pp_impr::Worm, partition::Worm, std::tuple<opBulla, op, bilinearPP>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<hedin_pp_impr::Worm, 1>, Mode, Value>(.5));
            markovChain.add(make_update<Reconnect<hedin_pp_impr::Worm, 2>, Mode, Value>(.5));
            
        }
        
        
        if(jParams.is(hedin_pp_imprsum::Worm::name())) {
            
            expansion::setup_updates<hedin_pp_imprsum::Worm, Mode>(jParams, data, state, markovChain);

            markovChain.add(make_update<InsertOps<partition::Worm, std::tuple<opBullaSum, op, bilinearPP>, hedin_pp_imprsum::Worm, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data),
                            make_update<RemoveOps<hedin_pp_imprsum::Worm, partition::Worm, std::tuple<opBullaSum, op, bilinearPP>, ut::sequence<0, 1, 2>>, Mode, Value>(1., jsx::empty_t(), data));

            markovChain.add(make_update<Reconnect<hedin_pp_imprsum::Worm, 1>, Mode, Value>(.5));
            markovChain.add(make_update<Reconnect<hedin_pp_imprsum::Worm, 2>, Mode, Value>(.5));
            
        }
    
       
        markovChain.finalize(state);
        
        
        mpi::cout << "End setting updates" << std::endl;
    };
    
}


#endif
