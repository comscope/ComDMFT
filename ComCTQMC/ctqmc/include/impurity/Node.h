#ifndef CTQMC_INCLUDE_IMPURITY_NODE_H
#define CTQMC_INCLUDE_IMPURITY_NODE_H

#include <iostream>
#include <random>

#include "Diagonal.h"
#include "Operators.h"
#include "../Utilities.h"

namespace imp {
    
    //Scheisse das muess verbessert werde chunt ja chei sau druus da ...

    template<typename NodeValue>
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

        int touched; NodeValue value;
    };
    
    
    template<typename NodeValue, typename F = void>
    struct Access {
        typedef Node<typename std::remove_cv<NodeValue>::type> NodeType;
        
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
        
        NodeValue* operator->() const { return &node_->value;};
        
    private:
        NodeType* node_;
        
        friend inline bool operator==(Access<NodeValue, F> const& lhs, Access<NodeValue, F> const& rhs) { return lhs.node_ == rhs.node_;};
        friend inline bool operator!=(Access<NodeValue, F> const& lhs, Access<NodeValue, F> const& rhs) { return lhs.node_ != rhs.node_;};
        
        friend F;
    };
    
    
    template<typename Mode, typename Value>
    struct NodeValue {
        NodeValue(Access<NodeValue const> node, itf::Operator<Value> const* op0, int flavor, EigenValues<Mode> const& eig) :
        op0(&get<Mode, Value>(*op0)), flavor(flavor),
        eig_(eig), node_(node),
        prop_(nullptr), propTry_(nullptr),
        ops_(new Operator<Mode, Value>*[2*node_.height()]), opsTry_(ops_ + node_.height()) {
            for(int l = 0; l < 2*node_.height(); ++l) ops_[l] = nullptr;
        };
        ~NodeValue() {
            for(int l = 0; l < 2*node_.height(); ++l) delete ops_[l]; 
            delete[] ops_; delete propTry_; delete prop_;
        };
        
        Operator<Mode, Value> const* const op0;
        int const flavor; 
        
        
        Propagator<Mode> const* prop() const {
            auto prop = node_.touched() ? prop_ : propTry_;
            if(prop == nullptr) throw std::runtime_error("imp::Value::prop: null pointer");
            return prop;
        };
        Operator<Mode, Value> const* op(int l) const {
            auto op = l < node_.touched() ? ops_[l] : opsTry_[l];
            if(op == nullptr) throw std::runtime_error("imp::Value::op: null pointer");
            return op;
        };
        std::unique_ptr<Operator<Mode, Value>> get_op(int l) const {
            auto& op = l < node_.touched() ? ops_[l] : opsTry_[l];
            if(op == nullptr) throw std::runtime_error("imp::Value::get_op: null pointer");
            auto temp = op; op = nullptr; return std::unique_ptr<Operator<Mode, Value>>(temp);
        };

        Propagator<Mode>* prop() {
            auto& prop = node_.touched() ? prop_ : propTry_;
            return prop ? prop : prop = new Propagator<Mode>(-(node_.next(0).key() - node_.key())*ut::beta()/ut::KeyMax, eig_);  //promotion stuff ...
        };
        Operator<Mode, Value>* op(int l) {
            auto& op = l < node_.touched() ? ops_[l] : opsTry_[l];
            return op ? op : op = new Operator<Mode, Value>(eig_);
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
        EigenValues<Mode> const& eig_;
        Access<NodeValue const> node_;
        
        Propagator<Mode>* prop_; Propagator<Mode>* propTry_;
        Operator<Mode, Value>** const ops_; Operator<Mode, Value>** const opsTry_;
    };
}

#endif


