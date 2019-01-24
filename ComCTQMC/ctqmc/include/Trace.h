#ifndef TRACE
#define TRACE

#include <iostream>
#include <random>

#include "Utilities.h"
#include "Updates.h"
#include "TraceElement.h"

// Achtung: nach construction 1) surviving, 2) prepare, 3) accept !!!!!!

namespace tr {
    
    enum class Flag { Pending, Accept, Reject };
	
	struct Trace {
        Trace() = delete;
        Trace(jsx::value const& jParams, EigenValues const& eig, Operators const& ops) :
        eig_(eig),
        ops_(ops), 
        urng_(std::mt19937(234), std::uniform_real_distribution<double>(.0, 1.)),
		prob_(jParams.is("skip-list probability") ? jParams("skip-list probability").real64() : .5),
		baseProb_(std::pow(prob_, (jParams.is("skip-list shift") ? jParams("skip-list shift").int64() : 0) + 1)),
		beta_(jParams("beta").real64()),
		maxHeight_(50),	
		height_(1), heightBackup_(1),
		size_(0), sizeBackup_(0),
        perm_(0), permBackup_(0),
        first_(new Element(-1, const_cast<Operator*>(&ops.identity()), maxHeight_, eig_)), //!!!!!!!!!!!!!! Porca miseria weg mit dem scheiss const_cast arsch vo chue.
		last_(nullptr),
        dimNorm_(eig.sectorNumber() + 1),
        it_(nullptr) {
			for(int l = 0; l < maxHeight_; ++l) first_->next(l) = last_;
            first_->distance(0) = Distance(ut::KeyMax, size_ + 1);
            inserted_.push_back(first_);  //!!!!!!!!!!! gits da chei besseri loesig ????? c.f. obe ..
            
            for(int s = eig.sectorNumber(); s; --s) dimNorm_[s] = std::log(eig_.at(s).dim());
			
            all_.begin = new SectorNorm*[(maxHeight_ + 1)*eig_.sectorNumber()];
            errorBounds_.begin = new ErrorBound[eig.sectorNumber() + 1];
            
            densityMatrix_ = std::unique_ptr<DensityMatrix>(new DensityMatrix(all_, eig_));
            
            int maxDim = 0; for(int s = eig.sectorNumber(); s; --s) maxDim = std::max(eig_.at(s).dim(), maxDim);
            bufferA_ = std::unique_ptr<Matrix>(new Matrix(maxDim*maxDim));
            bufferB_ = std::unique_ptr<Matrix>(new Matrix(maxDim*maxDim));
		}; 
        Trace(Trace const&) = delete;
        Trace(Trace&&) = delete;
        Trace& operator=(Trace const&) = delete;
        Trace& operator=(Trace&&) = delete;
        
		bool insert(ut::KeyType key, int flavor) {
            std::vector<Element*> elem(maxHeight_); Element* e = first_;
            std::vector<Distance> pos(maxHeight_); Distance p(0, 0);

			for(int l = height_; l--;) {
				for(; key > p.key + e->distance(l).key; e = e->next(l)) p += e->distance(l);
				elem[l] = e; pos[l] = p;
				++e->distance(l).entries;
			}

			if(key == p.key + e->distance(0).key) {
				for(int l = height_; l--;) --elem[l]->distance(l).entries;
                return false;
			} 
            ++p.entries; p.key = key; ++size_;
			
			int const newHeight = randomHeight();
			for(int l = height_; l < newHeight; ++l) { 
				elem[l] = first_; pos[l] = Distance(0, 0);
				first_->distance(l) = Distance(ut::KeyMax, size_ + 1);  
			}

			e = new Element(flavor, const_cast<Operator*>(&ops_.at(flavor)), newHeight, eig_);
			inserted_.push_back(e);
			
			for(int l = 0; l < newHeight; ++l) {
				e->next(l) = elem[l]->next(l);			
				elem[l]->next(l) = e;
				
				Distance const diff = p - pos[l];
				e->distance(l) = elem[l]->distance(l) - diff;
				elem[l]->distance(l) = diff;
			}

			if(newHeight > height_) height_ = newHeight; 
			
			elem[0]->touch(touched_, 0);
			for(int l = 1; l < height_; ++l)
				if(elem[l - 1] != elem[l]) elem[l]->touch(touched_, l);
            
            perm_ += p.entries;
            return true;
		};
		
		void erase(ut::KeyType key) {
            std::vector<Element*> elem(maxHeight_); Element* e = first_;
			Distance p(0, 0);
			
			for(int l = height_; l--;) {
				for(; key > p.key + e->distance(l).key; e = e->next(l)) p += e->distance(l);
				elem[l] = e; --e->distance(l).entries;
			}
			
			if(key != p.key + e->distance(0).key) {
				for(int l = height_; l--;) ++elem[l]->distance(l).entries;
                throw std::runtime_error("Tr: key not found !");
			}
			
			e = e->next(0);
			erased_.push_back(e);
			++p.entries; --size_;
			
			for(int l = 0; l < e->height(); ++l) {
				elem[l]->next(l) = e->next(l);
				elem[l]->distance(l) += e->distance(l);
			};
			
			elem[0]->touch(touched_, 0);
			for(int l = 1; l < height_; ++l)
				if(elem[l - 1] != elem[l]) elem[l]->touch(touched_, l);
			
			for(; height_ > 1 && first_->next(height_ - 1) == last_; --height_)
                ;
			
            perm_ += p.entries;  //modulo hier ziemlich sicher nicht notwending ... guck standard wegen wraparound .... pass auf in sign() falls nicht notwending ...
		};
		
        Flag surviving() {
            densityMatrixTry_ = std::unique_ptr<DensityMatrix>(new DensityMatrix(all_, eig_));
            
            if(size_)
                map(first_, height_, all_);
            else
                first_->prop()->add(all_);
            
            if(densityMatrixTry_->surviving(errorBounds_)) return Flag::Pending;
            
            return Flag::Reject;
        };
        
        void prepare(za::Zahl x) {
            traceThresh_ = za::abs(x*densityMatrix_->Z());
            
            for(ErrorBounds::iterator it = errorBounds_.begin; it != errorBounds_.end; ++it) it->norm += dimNorm_[it->sector];
            
            std::sort(errorBounds_.begin, errorBounds_.end, &compare);
            
            for(ErrorBounds::iterator it = errorBounds_.begin; it != errorBounds_.end; ++it) it->bound = za::pow(it->norm);
            for(ErrorBounds::iterator it = errorBounds_.end - 1; it != errorBounds_.begin; --it) (it - 1)->bound += it->bound;
            errorBounds_.end->bound = za::Zahl(.0);
            
            it_ = errorBounds_.begin;
        }
        
        Flag decide() {
            if(za::abs(densityMatrixTry_->Z()) + it_->bound <= traceThresh_)
                return Flag::Reject;
            else if(it_->bound <= std::numeric_limits<double>::epsilon()*za::abs(densityMatrixTry_->Z()))
                return Flag::Accept;
            
            if(size_)
                prod(densityMatrixTry_.get(), first_, height_, it_->sector, bufferA_.get(), bufferB_.get());
            else
                densityMatrixTry_->mat(it_->sector, -beta_, eig_.at(it_->sector));
            
            trace(&densityMatrixTry_->Z(it_->sector), &densityMatrixTry_->Z(), densityMatrixTry_->mat(it_->sector));
            
            ++it_; return Flag::Pending;
        };
		
		
		void accept() {
			heightBackup_ = height_;
			
            for(auto ptr : touched_) ptr->accept();
            for(auto ptr : inserted_) ptr->accept();
            for(auto ptr : erased_) delete ptr;
			
			touched_.clear(); inserted_.clear(); erased_.clear();
						
            densityMatrix_ = std::move(densityMatrixTry_); 
			
			sizeBackup_ = size_; permBackup_ = perm_;
		};
		
		void reject() {
			height_ = heightBackup_;
			
            for(auto ptr : touched_) ptr->reject();
            for(auto ptr : inserted_) delete ptr;
			
			touched_.clear(); inserted_.clear(); erased_.clear();
			
			densityMatrixTry_.reset();
            
            size_ = sizeBackup_; perm_ = permBackup_;
		}; 
											 
		int sign() const {
            return (perm_%2 ? -1 : 1)*(densityMatrix_->Z() <= .0 ? -1 : 1);
        };
        int size() const {
            return size_;
        };
		DensityMatrix const& densityMatrix() const {
            return *densityMatrix_;
        };
		
		~Trace() {
            reject();
            
            delete[] errorBounds_.begin;
            delete[] all_.begin;

			Element* it = first_;		
			while(it != last_) {
				Element* e = it;
				it = it->next(0);
				delete e;
			}
		};
    private:
        EigenValues const& eig_;
        Operators const& ops_;

		ut::RandomNumberGenerator<std::mt19937, std::uniform_real_distribution<double> > urng_;
		
		double const prob_; 
		double const baseProb_;
		double const beta_;
		int const maxHeight_;	
		int height_, heightBackup_;
		int size_, sizeBackup_;
        unsigned int perm_, permBackup_;

		Element* const first_; Element* const last_;		

		std::vector<Element*> touched_;
		std::vector<Element*> inserted_;
		std::vector<Element*> erased_;
        
        std::vector<double> dimNorm_;
		
        SectorNormPtrs all_; ErrorBounds errorBounds_;
        
        ErrorBounds::iterator it_; za::Zahl traceThresh_;
        
        std::unique_ptr<DensityMatrix> densityMatrix_;
        std::unique_ptr<DensityMatrix> densityMatrixTry_;
		
        std::unique_ptr<Matrix> bufferA_;
        std::unique_ptr<Matrix> bufferB_;
		
		int randomHeight() {
			int l = 1; double u = urng_();
			for(double p = baseProb_; u < p && l < maxHeight_ - 1; ++l, p *= prob_); 
			return l;
		};
		
		Distance distance(Element* elem, int level) {
			Element* const end = elem->next(level); Distance dist(0, 0);
			for(; elem != end; elem = elem->next(level - 1)) dist += elem->distance(level - 1);
			return dist;
		};

		static void map(Element* elem, int level, SectorNormPtrs& requested) {
			Element* const last = elem->next(level);
            while(elem->next(level) == last) --level;
			
			do {
				if(elem->next(level) != elem->next(0)) {
					SectorNormPtrs missing;
					if(elem->op(level)->missing(missing, requested)) map(elem, level, missing);
					elem->op(level)->map(requested);
				} else {
					elem->op(0)->map(requested);
					elem->prop()->add(requested); 
				} 
				elem = elem->next(level);
			} while(elem != last);
		}
        
        //!!!!!!!!!!!!!!!!!!!! Das gaht en stuck chuerzer da ......
        static void prod(BlocMatrixBase* const dest, Element* elem, int level, int sec, Matrix* A, Matrix* B) {
            Element* const last = elem->next(level);
            while(elem->next(level) == last) --level;
            
            Operator* op; int s = sec;
            for(Element* e = elem; e != last; e = e->next(level)) {
                if(e->next(level) != e->next(0)) {
                    op = e->op(level);
                    if(!op->is(s)) {
                        prod(op, e, level, s, A, B); norm(&op->norm(s), op->mat(s));
                    }
                } else
                    op = e->op(0);
                s = op->sec(s);
            }
            
            Matrix* previous; s = sec;
            if(elem->next(level) != elem->next(0)) {
                op = elem->op(level); previous = &op->mat(s); s = op->sec(s);
            } else {
                op = elem->op(0); previous = &op->mat(s); s = op->sec(s);
                copyEvolveL(elem->prop()->at(s), *A, *previous);
                previous = A;
            }
            elem = elem->next(level);
            
            while(elem->next(level) != last) {
                if(elem->next(level) != elem->next(0)) {
                    op = elem->op(level); mult(*B, op->mat(s), *previous); s = op->sec(s);
                } else {
                    op = elem->op(0); mult(*B, op->mat(s), *previous); s = op->sec(s);
                    evolveL(elem->prop()->at(s), *B);
                }
                elem = elem->next(level); std::swap(A, B); previous = A;
            }
            
            if(elem->next(level) != elem->next(0)) {
                op = elem->op(level); mult(dest->mat(sec, op->mat(s).I()*previous->J()), op->mat(s), *previous);
            } else {
                op = elem->op(0); mult(dest->mat(sec, op->mat(s).I()*previous->J()), op->mat(s), *previous);
                evolveL(elem->prop()->at(op->sec(s)), dest->mat(sec));
            }
        };
	};
    
    template<typename> struct Updates {};
    
    template<>
    struct Updates<up::InsertTwo> {
        Flag surviving(up::InsertTwo const& u, Trace& trace) {
            if(!trace.insert(u.keyR(), u.flavorR)) return Flag::Reject;
            if(!trace.insert(u.keyL(), u.flavorL)) return Flag::Reject;
            return trace.surviving();
        };
        void prepare(Trace& trace, za::Zahl x) { trace.prepare(x);};
        Flag decide(Trace& trace) { return trace.decide();};
        void accept(Trace& trace) { trace.accept();};
        void reject(Trace& trace) { trace.reject();};
    };
    
    template<>
    struct Updates<up::EraseTwo> {
        Flag surviving(up::EraseTwo const& u, Trace& trace) {
            trace.erase(u.keyL());
            trace.erase(u.keyR());
            return trace.surviving();
        };
        void prepare(Trace& trace, za::Zahl x) { trace.prepare(x);};
        Flag decide(Trace& trace) { return trace.decide();};
        void accept(Trace& trace) { trace.accept();};
        void reject(Trace& trace) { trace.reject();};
    };
}
	
#endif 


