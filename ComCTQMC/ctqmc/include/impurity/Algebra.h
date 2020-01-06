#ifndef IMPURITY_ALGEBRA_H
#define IMPURITY_ALGEBRA_H

#include <cmath>
#include <iostream>
#include <vector>

#include "../Utilities.h"
#include "../../../include/JsonX.h"
#include "../../../include/io/Matrix.h"


namespace imp {
    
    namespace itf {
        
        struct Batcher {
            virtual int is_ready() = 0;
            virtual void launch() = 0;
            virtual ~Batcher() = default;
        };
        
    };
    
    template<typename Alloc> struct Batcher;
    
    template<typename Alloc> Batcher<Alloc>& get(itf::Batcher&);
    
    
    template<typename Alloc>
	struct Energies {
        using data_type = typename Alloc::template Data<double>; //An Haesslichkeit schwer z'uebertreffe ......
        
        Energies() = delete;
        Energies(jsx::value const& jParams, std::vector<double> const& eig);
        Energies(Energies const&) = delete;
        Energies(Energies&&) = delete;
        Energies& operator=(Energies const&) = delete;
        Energies& operator=(Energies&&) = delete;
        ~Energies();
        
		int const& dim0() const { return dim0_;};
        int const& dim() const { return dim_;}
        double const& ln_dim() const { return ln_dim_;};
        data_type& data() { return data_;}
        data_type const& data() const { return data_;}
		double const& min() const { return min_;}
        
	private:
        int const dim0_;
		int const dim_;
        double const ln_dim_;
		data_type data_;
		double min_;
    };

    template<typename Alloc>
	struct Vector {
        using data_type = typename Alloc::template Data<double>;
        
        Vector() = delete;
		Vector(double time, Energies<Alloc> const& energies);
        Vector(Vector const&) = delete;
        Vector(Vector&&) = delete;
        Vector& operator=(Vector const&) = delete;
        Vector& operator=(Vector&&) = delete;
        ~Vector();
        
        double const& time() const { return time_;};
		double const& exponent() const { return exponent_;}
        Energies<Alloc> const& energies() const { return energies_;};  
        data_type& data() { return data_;}
		data_type const& data() const { return data_;}
		
	private:
        double const time_;
		double const exponent_;
        Energies<Alloc> const& energies_;    /// This is a bit shity, eig_ (data_) is not used in host (device)
		data_type data_;
	};
	
    template<typename Alloc>
	struct Matrix {
        using data_type = typename Alloc::template Data<double>;
        
        struct Identity { Identity(int d) : dim(d) {}; int const dim;};
        struct Zero { Zero(int d) : dim(d) {}; int const dim;};
        
        Matrix() = delete;
        Matrix(int size);
        Matrix(Identity const& identity);
        Matrix(Zero const& zero);
        Matrix(int I, int J, io::rmat const& mat);
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = delete;
        Matrix& operator=(Matrix const&) = delete;
        Matrix& operator=(Matrix&&) = delete;
        ~Matrix();
        
        int& I() { return I_;}
        int& J() { return J_;}
		int const& I() const { return I_;}
		int const& J() const { return J_;}
        data_type& data() { return data_;}
        data_type const& data() const { return data_;}
        double& exponent() { return exponent_;}
        double const& exponent() const { return exponent_;}
    
	protected:        
        int I_, J_;
        data_type data_;
		double exponent_;
	};
}


#endif
