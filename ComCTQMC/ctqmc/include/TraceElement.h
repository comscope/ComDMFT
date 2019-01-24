#ifndef TRACEELEMENT
#define TRACEELEMENT

#include <iostream>
#include <algorithm>
#include <vector>
#include "Utilities.h"
#include "TraceUtilities.h"

///////////////////////////////////////////////
//  skip list mit unique_ptr's !!!!
///////////////////////////////////////////////

namespace tr {

	struct Distance {
		Distance() {};
		Distance(ut::KeyType k, int e) : key(k), entries(e) {};
		Distance& operator+=(Distance const& d) {
			key += d.key;
			entries += d.entries;
			return *this;
		};		
		ut::KeyType key; int entries;
	};
	
	Distance operator-(Distance const& lhs, Distance const& rhs) {
		return Distance(lhs.key - rhs.key, lhs.entries - rhs.entries);
	}
	
	struct Element {
        Element() = delete;
		Element(int flavor, Operator* op, int height, EigenValues const& eig) :
        eig_(eig),
		flavor_(flavor),
		height_(height),
		touchedLevel_(0),
		next_(new Element*[2*height]), 
		nextBackup_(next_ + height),
		distance_(new Distance[2*height]),
		distanceBackup_(distance_ + height),
		prop_(0),
		propTry_(0),
		ops_(new Operator*[2*height]),
		opsTry_(ops_ + height) { 
			for(int l = 0; l < 2*height_; ++l) ops_[l] = 0; 
			ops_[0] = opsTry_[0] = op;
		};
        Element(Element const&) = delete;
        Element(Element&&) = delete;
        Element& operator=(Element const&) = delete;
        Element& operator=(Element&&) = delete;
 		
		int flavor() const { return flavor_;};
		int height() const { return height_;};
		Element*& next(int l) { return next_[l];};
		Distance& distance(int l) { return distance_[l];};
		
		Propagator* prop() { 
			Propagator*& prop = touchedLevel_ ? prop_ : propTry_;
			return prop ? prop : prop = new Propagator(-distance_[0].key*ut::beta()/ut::KeyMax, eig_);
		};
		
		Operator* op(int l) { 	
			Operator*& op = l < touchedLevel_ ? ops_[l] : opsTry_[l];
			return op ? op : op = new Operator(eig_);
		};                     
		
		void touch(std::vector<Element*>& touched, int l) {
			if(touchedLevel_ == height_) touched.push_back(this);
			touchedLevel_ = std::min(l, touchedLevel_);	
		};
		
		void accept() {
			std::memcpy(nextBackup_ + touchedLevel_, next_ + touchedLevel_, (height_ - touchedLevel_)*sizeof(Element*));
			std::memcpy(distanceBackup_ + touchedLevel_, distance_ + touchedLevel_, (height_ - touchedLevel_)*sizeof(Distance));
			
			if(!touchedLevel_) { delete prop_; prop_ = propTry_; propTry_ = 0;}
			for(int l = std::max(touchedLevel_, 1); l < height_; ++l) { delete ops_[l]; ops_[l] = opsTry_[l]; opsTry_[l] = 0;}
			touchedLevel_ = height_;
		};
		
		void reject() {
			std::memcpy(next_ + touchedLevel_, nextBackup_ + touchedLevel_, (height_ - touchedLevel_)*sizeof(Element*));
			std::memcpy(distance_ + touchedLevel_, distanceBackup_ + touchedLevel_, (height_ - touchedLevel_)*sizeof(Distance));
			
			if(!touchedLevel_) { delete propTry_; propTry_ = 0;}
			for(int l = std::max(touchedLevel_, 1); l < height_; ++l) { delete opsTry_[l]; opsTry_[l] = 0;}
			touchedLevel_ = height_;
		};
		
		~Element() {
			ops_[0] = opsTry_[0] = 0;
			for(int l = 0; l < 2*height_; ++l)  delete ops_[l];
			delete[] ops_;
			
			delete propTry_;
			delete prop_; 
			
			delete[] distance_;
			delete[] next_;
		};
	private:
        EigenValues const& eig_;
        
		int const flavor_;
		int const height_;	
		
		int touchedLevel_;
		Element** const next_; Element** const nextBackup_;
		Distance* const distance_; Distance* const distanceBackup_;
		Propagator* prop_; Propagator* propTry_; 
		Operator** const ops_; Operator** const opsTry_;
	};
}
	
#endif 