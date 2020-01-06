#ifndef IMPURITY_NODE_H
#define IMPURITY_NODE_H

#include <iostream>
#include <random>

#include "Diagonal.h"
#include "Operators.h"
#include "../Utilities.h"

namespace imp {

    template<typename ValueType>
    struct Node {
        template<typename... Args>
        Node(ut::KeyType key, int height, Args&&... args) :
        key(key), height(height),
        next(new Node*[2*height]),
        entries(new int[2*height]),
        touched(0),
        value(this, std::forward<Args>(args)...) {
        };
        ~Node() {
            delete[] entries; delete[] next;
        };
        
        ut::KeyType const key; int const height;
        Node** const next; int* const entries;
        
        void accept() {
            std::memcpy(next + height + touched, next + touched, (height - touched)*sizeof(Node*));
            std::memcpy(entries + height + touched, entries + touched, (height - touched)*sizeof(int));
            value.accept(); touched = height;
        };
        void reject() {
            std::memcpy(next + touched, next + height + touched, (height - touched)*sizeof(Node*));
            std::memcpy(entries + touched, entries + height + touched, (height - touched)*sizeof(int));
            value.reject(); touched = height;
        };
        bool touch(int level) {
            bool temp = (touched == height);
            touched = std::min(level, touched);
            return temp;
        };

        int touched; ValueType value;
    };
    
    
    template<typename ValueType, typename F = void>
    struct Access {
        typedef Node<typename std::remove_cv<ValueType>::type> NodeType;
        
        Access() = default;
        Access(NodeType* node) : node_(node) {};
        Access(Access const&) = default;
        Access(Access&&) = default;
        Access& operator=(Access const&) = default;
        Access& operator=(Access&&) = default;
        ~Access() = default;
        
        int height() const { return node_->height;};
        ut::KeyType key() const { return node_->key;};
        Access next(int level) const { return node_->next[level];};
        int entries(int level) const { return node_->entries[level];};
        int touched() const { return node_->touched;};
        
        ValueType* operator->() const { return &node_->value;};
        
    private:
        NodeType* node_;
        
        friend inline bool operator==(Access<ValueType, F> const& lhs, Access<ValueType, F> const& rhs) { return lhs.node_ == rhs.node_;};
        friend inline bool operator!=(Access<ValueType, F> const& lhs, Access<ValueType, F> const& rhs) { return lhs.node_ != rhs.node_;};
        
        friend F;
    };
    
    
    template<typename Alloc>
    struct Value {
        Value(Access<Value const> node, Operator<Alloc> const* op0, int flavor, double* ptr, EigenValues<Alloc> const& eig) :
        op0(op0), flavor(flavor), ptr(ptr),
        eig_(eig), node_(node),
        prop_(nullptr), propTry_(nullptr),
        ops_(new Operator<Alloc>*[2*node_.height()]), opsTry_(ops_ + node_.height()) {
            for(int l = 0; l < 2*node_.height(); ++l) ops_[l] = nullptr;
        };
        ~Value() {
            for(int l = 0; l < 2*node_.height(); ++l) delete ops_[l]; delete[] ops_;
            delete propTry_; delete prop_;
        };
        
        Operator<Alloc> const* const op0;
        int const flavor; double* const ptr;
        
        
        Propagator<Alloc> const* prop() const {
            auto prop = node_.touched() ? prop_ : propTry_;
            if(prop == nullptr) throw std::runtime_error("imp::Value::prop: null pointer");
            return prop;
        };
        Operator<Alloc> const* op(int l) const {
            auto op = l < node_.touched() ? ops_[l] : opsTry_[l];
            if(op == nullptr) throw std::runtime_error("imp::Value::op: null pointer");
            return op;
        };
        std::unique_ptr<Operator<Alloc>> get_op(int l) const {
            auto& op = l < node_.touched() ? ops_[l] : opsTry_[l];
            if(op == nullptr) throw std::runtime_error("imp::Value::get_op: null pointer");
            auto temp = op; op = nullptr; return std::unique_ptr<Operator<Alloc>>(temp);
        };

        Propagator<Alloc>* prop() {
            auto& prop = node_.touched() ? prop_ : propTry_;
            return prop ? prop : prop = new Propagator<Alloc>(-(node_.next(0).key() - node_.key())*ut::beta()/ut::KeyMax, eig_);  //promotion stuff ...
        };
        Operator<Alloc>* op(int l) {
            auto& op = l < node_.touched() ? ops_[l] : opsTry_[l];
            return op ? op : op = new Operator<Alloc>(eig_);
        };
        
        void accept() {
            if(!node_.touched()) {
                delete prop_; prop_ = propTry_; propTry_ = nullptr;
            }
            for(int l = std::max(node_.touched(), 1); l < node_.height(); ++l) {
                delete ops_[l]; ops_[l] = opsTry_[l]; opsTry_[l] = nullptr;
            }
        };
        void reject() {
            if(!node_.touched()) {
                delete propTry_; propTry_ = nullptr;
            }
            for(int l = std::max(node_.touched(), 1); l < node_.height(); ++l) {
                delete opsTry_[l]; opsTry_[l] = nullptr;
            }
        };
        
    private:
        EigenValues<Alloc> const& eig_;
        Access<Value const> node_;
        
        Propagator<Alloc>* prop_; Propagator<Alloc>* propTry_;
        Operator<Alloc>** const ops_; Operator<Alloc>** const opsTry_;
    };
}

#endif


