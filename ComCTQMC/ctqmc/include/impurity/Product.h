#ifndef CTQMC_INCLUDE_IMPURITY_PRODUCT_H
#define CTQMC_INCLUDE_IMPURITY_PRODUCT_H

#include <iostream>
#include <random>

#include "Algebra.h"
#include "Node.h"
#include "../Utilities.h"
#include "../../../include/JsonX.h"

namespace imp {
    
    namespace itf {
        
        template<typename Value>
        struct Product {
            virtual int size() const = 0;
            virtual int height() const = 0;
            virtual bool insert(ut::KeyType const key, int flavor) = 0;
            virtual void erase(ut::KeyType key) = 0;
            virtual int accept() = 0;
            virtual void reject() = 0;
            virtual ~Product() = default;
        };
        
    };
    
    
    // Todo: remove flavor from node 
    
    template<typename Mode, typename Value>
    struct Product : itf::Product<Value> {
        using NodeType  = Node<NodeValue<Mode, Value>>;
        
        using AccessType  = Access<NodeValue<Mode, Value>, Product>;
        using CAccessType = Access<NodeValue<Mode, Value> const, Product>;
        
        Product() = delete;
        Product(jsx::value const& jParams, itf::EigenValues const& eig, itf::Operator<Value> const& ide, itf::Operators<Value> const& ops) :
        eig_(get<Mode>(eig)), ide_(get<Mode>(ide)), ops_(get<Mode>(ops)),
        urng_(std::mt19937(234), std::uniform_real_distribution<double>(.0, 1.)),
		prob_(jParams.is("skip-list probability") ? jParams("skip-list probability").real64() : .5),
		baseProb_(std::pow(prob_, (jParams.is("skip-list shift") ? jParams("skip-list shift").int64() : 0) + 1)),
		maxHeight_(50),
		height_(1), heightBackup_(1),
		size_(0), sizeBackup_(0),
        sign_(1),
        first_(new NodeType(         0, maxHeight_ + 1,   &ide_, -1, eig_)),
        last_( new NodeType(ut::KeyMax,              0, nullptr, -1, eig_)) {
            for(int l = 0; l <= height_; ++l) { first_->next[l] = last_; first_->entries[l] = size_ + 1;}
            first_->accept();
            
            all_.begin = new SectorNorm*[(maxHeight_ + 2)*eig.sectorNumber()];
            
            int maxDim = 0; for(int s = eig.sectorNumber(); s; --s) maxDim = std::max(eig_.at(s).dim(), maxDim);
            bufferA_.reset(new Matrix<Mode, Value>(maxDim*maxDim));
            bufferB_.reset(new Matrix<Mode, Value>(maxDim*maxDim));
            bufferC_.reset(new Matrix<Mode, Value>(maxDim*maxDim));
		};
        Product(Product const&) = delete;
        Product(Product&&) = delete;
        Product& operator=(Product const&) = delete;
        Product& operator=(Product&&) = delete;
        ~Product() {
            reject();
            
            delete[] all_.begin;

            auto it = first_;
            while(it != last_) {
                auto temp = it;
                it = it->next[0];
                delete temp;
            }
            delete last_;
        };

        int size() const { return size_;};
        int height() const { return height_;};

        CAccessType first() const { return first_;};
        CAccessType last() const { return last_;};
        
        bool insert(ut::KeyType const key, int flavor) {
            return last() != insert_impl(key, random_height(), &ops_.at(flavor), flavor, eig_);
        };
        // Todo: remove flavor entry
        CAccessType insert(ut::KeyType const key, itf::Operator<Value> const* op, int flavor) {
            return insert_impl(key, random_height(), &get<Mode>(*op), flavor, eig_);
        };
        CAccessType insert(ut::KeyType const key, itf::Operator<Value> const* op, int flavor, int height) {
            if(height > maxHeight_ + 1) throw std::runtime_error("imp::Product::insert: invalid height");
            return insert_impl(key, height, &get<Mode>(*op), flavor, eig_);
        };
        
        void erase(ut::KeyType key) {
            erase_impl(key);
        };
        
        std::vector<int> map(CAccessType cbegin, int level, std::vector<int> const& sectors) {
            auto begin = AccessType(cbegin.node_); std::vector<int> result;

            begin->op(level)->assign(sectors, all_);
            begin.next(level) != begin.next(0) ? map_impl(begin, level, all_) : begin->prop()->add(all_);
            for(auto sec : sectors)
                result.push_back(cbegin->op(level)->map(sec).sector);

            return result;
        };
        void multiply(CAccessType cbegin, int level, int sec, itf::Batcher<Value>& batcher) {
            auto begin = AccessType(cbegin.node_); if(begin->op(level)->isMat(sec)) return;

            if(begin.next(level) != begin.next(0)) {
                if(size_ + 2 > ptrMat_.size()) { ptrMat_.resize(size_ + 2); ptrProp_.resize(size_ + 2);}
                multiply_impl(begin, level, sec,  ptrMat_.data(), ptrProp_.data(), bufferA_.get(), bufferB_.get(), batcher);
            } else {
                copyEvolveL(begin->op(level)->mat(sec, ide_.mat(sec).I()*ide_.mat(sec).J()), begin->prop()->at(sec), begin->op0->mat(sec), batcher);
                auto& map = begin->op(level)->set_map(sec); map.sector = begin->op0->map(sec).sector; norm(&map.norm, begin->op(level)->mat(sec), batcher);
            }
        };
		
		int accept() {
			heightBackup_ = height_;
			
            for(auto ptr : touched_) ptr->accept();
            for(auto ptr : inserted_) ptr->accept();
            for(auto ptr : erased_) delete ptr;
			
			touched_.clear(); inserted_.clear(); erased_.clear();

            sizeBackup_ = size_;
            
            int temp = sign_; sign_ = 1; return temp;
        };
		void reject() {
			height_ = heightBackup_;
			
            for(auto ptr : touched_) ptr->reject();
            for(auto ptr : inserted_) delete ptr;
			
			touched_.clear(); inserted_.clear(); erased_.clear();

            size_ = sizeBackup_;
            
            sign_ = 1;
		};
        
        Matrix<Mode, Value>& bufferA() { return *bufferA_;};
        Matrix<Mode, Value>& bufferB() { return *bufferB_;};
        Matrix<Mode, Value>& bufferC() { return *bufferC_;};

    private:
        EigenValues<Mode> const& eig_;
        Operator<Mode, Value> const& ide_;
        Operators<Mode, Value> const& ops_;
        
		ut::RandomNumberGenerator<std::mt19937, std::uniform_real_distribution<double> > urng_;
		
		double const prob_; 
		double const baseProb_;
		int const maxHeight_;	
		int height_, heightBackup_;
		int size_, sizeBackup_;
        int sign_;

		NodeType* const first_; NodeType* const last_;		

		std::vector<NodeType*> touched_, inserted_, erased_;
        
        SectorNormPtrs all_;
        
        std::vector<Matrix<Mode, Value> const*> ptrMat_;
        std::vector<Vector<Mode> const*> ptrProp_;
        std::unique_ptr<Matrix<Mode, Value>> bufferA_, bufferB_, bufferC_;


		int random_height() {
			int h = 1; double u = urng_();
			for(double p = baseProb_; u < p && h < maxHeight_; ++h, p *= prob_);
			return h;
        };
        
        template<typename... Args>
        NodeType* insert_impl(ut::KeyType const key, int const newHeight, Args&&... args) {
            if(!(0 <= key && key <= ut::KeyMax)) throw std::runtime_error("imp::Product::insert_impl: invalid key");

            std::vector<NodeType*> node(maxHeight_ + 1); NodeType* n = first_;
            std::vector<int> pos(maxHeight_ + 1); int p = 0;

            for(int l = std::min(height_, maxHeight_); l >= 0; --l) {
                for(; key > n->next[l]->key; n = n->next[l]) p += n->entries[l];
                node[l] = n; pos[l] = p; ++n->entries[l];
            }

            if(key == n->next[0]->key) {
                for(int l = std::min(height_, maxHeight_); l >= 0; --l) --node[l]->entries[l];
                return last_;
            }
            ++p; ++size_;

            for(int l = height_ + 1; l <= std::min(newHeight, maxHeight_); ++l) {
                node[l] = first_; pos[l] = 0; first_->next[l] = last_; first_->entries[l] = size_ + 1;
            }

            n = new NodeType(key, newHeight, std::forward<Args>(args)...); 
            inserted_.push_back(n);

            for(int l = 0; l < newHeight; ++l) {
                n->next[l] = node[l]->next[l];
                node[l]->next[l] = n;

                int const diff = p - pos[l];
                n->entries[l] = node[l]->entries[l] - diff;
                node[l]->entries[l] = diff;
            }

            if(newHeight > height_) height_ = newHeight;
            
            if(node[0]->touch(0)) touched_.push_back(node[0]);
            for(int l = 1; l <= std::min(height_, maxHeight_); ++l)
                if(node[l - 1] != node[l])
                    if(node[l]->touch(l)) touched_.push_back(node[l]);

            sign_ *= p%2 ? -1 : 1; return n;
        };
        
        void erase_impl(ut::KeyType const key) {
            if(!(0 < key && key < ut::KeyMax)) throw std::runtime_error("imp::Product::erase_impl: invalid key");
            
            std::vector<NodeType*> node(maxHeight_ + 1); NodeType* n = first_;
            int p = 0;
            
            for(int l = std::min(height_, maxHeight_); l >= 0; --l) {
                for(; key > n->next[l]->key; n = n->next[l]) p += n->entries[l];
                node[l] = n; --n->entries[l];
            }
            n = n->next[0];
            
            if(key != n->key) throw std::runtime_error("imp::Product::erase_impl: key not found !");
            
            erased_.push_back(n);
            ++p; --size_;
            
            for(int l = 0; l < n->height; ++l) {
                node[l]->next[l] = n->next[l];
                node[l]->entries[l] += n->entries[l];
            };
            
            if(node[0]->touch(0)) touched_.push_back(node[0]);
            for(int l = 1; l <= std::min(height_, maxHeight_); ++l)
                if(node[l - 1] != node[l])
                    if(node[l]->touch(l)) touched_.push_back(node[l]);

            for(; height_ > 1 && first_->next[height_ - 1] == last_; --height_)
                ;
            
            sign_ *= p%2 ? -1 : 1;  //modulo hier ziemlich sicher nicht notwending ... guck standard wegen wraparound .... pass auf in sign() falls nicht notwending ...
        };
        
        
        static void map_impl(AccessType it, int level, SectorNormPtrs& requested) {
            auto const last = it.next(level);
            while(it.next(level) == last) --level;

            do {
                if(it.next(level) != it.next(0)) {
                    SectorNormPtrs missing;
                    if(it->op(level)->missing(missing, requested)) map_impl(it, level, missing);
                    it->op(level)->map(requested);
                } else {
                    it->op0->map(requested);
                    it->prop()->add(requested);
                }
                it = it.next(level);
            } while(it != last);
        }
        
        static void multiply_impl(AccessType const begin, int const level, int const sec, Matrix<Mode, Value> const** const mat, Vector<Mode> const** const prop, Matrix<Mode, Value>* A, Matrix<Mode, Value>* B, itf::Batcher<Value>& batcher) {
            auto const end = begin.next(level);
            auto l = level; while(begin.next(l) == end) --l;

            auto s = sec; auto m = mat; auto p = prop;
            for(auto it = begin; it != end; it = it.next(l))
                if(it.next(l) != it.next(0)) {
                    if(!it->op(l)->isMat(s)) multiply_impl(it, l, s, m, p, A, B, batcher);
                    *m++ = &it->op(l)->mat(s); s = it->op(l)->map(s).sector; *p++ = nullptr;
                } else {
                    *m++ = &it->op0->mat(s); s = it->op0->map(s).sector; *p++ = &it->prop()->at(s);
                }

            *m = nullptr; m = mat; p = prop;
            if(*p != nullptr) {
                copyEvolveL(*A, **p, **m, batcher);
                *m = A; std::swap(A, B);
            }
            while(*(++m + 1) != nullptr) {
                mult(*A, **m, **(m - 1), batcher);
                if(*++p != nullptr) evolveL(**p, *A, batcher);
                *m = A; std::swap(A, B);
            }
            mult(begin->op(level)->mat(sec, (*m)->I()*(*(m - 1))->J()), **m, **(m - 1), batcher);
            if(*++p != nullptr) evolveL(**p, begin->op(level)->mat(sec), batcher);
            
            auto& map = begin->op(level)->set_map(sec);
            map.sector = s; norm(&map.norm, begin->op(level)->mat(sec), batcher);
        };
    };
    
    template<typename Mode, typename Value> Product<Mode, Value>& get(itf::Product<Value>& productItf) {
        return static_cast<Product<Mode, Value>&>(productItf);
    };
    
    template<typename Mode, typename Value> Product<Mode, Value> const& get(itf::Product<Value> const& productItf) {
        return static_cast<Product<Mode, Value> const&>(productItf);
    };
}

#endif


