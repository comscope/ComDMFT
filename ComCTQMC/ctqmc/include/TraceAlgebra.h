#ifndef TRACEALGEBRA
#define TRACEALGEBRA

#include <cmath>
#include <iostream>
#include <vector>

#include "Zahl.h"

#include "../../include/BlasLapack.h"
#include "../../include/JsonX.h"
#include "../../include/IO.h"


namespace tr {
	
	struct Energies {
        Energies() = delete;
        Energies(jsx::value const& jParams, std::vector<double> const& eig);
        Energies(Energies const&) = delete;
        Energies(Energies&&) = delete;
        Energies& operator=(Energies const&) = delete;
        Energies& operator=(Energies&&) = delete;
		int const& dim0() const { return dim0_;};
        int const& dim() const { return dim_;}
		data_ptr data() const { return data_;}
		double const& min() const { return min_;}		
		~Energies();
	private:
        int const dim0_;
		int const dim_;
		data_ptr data_;
		double min_;
	};

	struct Vector {
        Vector() = delete;
		Vector(double time, Energies const& eig);
        Vector(Vector const&) = delete;
        Vector(Vector&&) = delete;
        Vector& operator=(Vector const&) = delete;
        Vector& operator=(Vector&&) = delete;
        double const& time() const { return time_;};
		double const& exponent() const { return exponent_;}
		data_ptr data() const { return data_;}
		~Vector();
	private:
        //int const dim_;
        double const time_;
		double const exponent_;
		data_ptr data_;
	};
	
	struct Matrix {
        struct Identity { Identity(int d) : dim(d) {}; int const dim;};
        struct Zero { Zero(int d) : dim(d) {}; int const dim;};
        
        Matrix() = delete;
        Matrix(int size);
        Matrix(Identity const& identity);
        Matrix(int I, int J, io::rmat const& mat);
        Matrix(Zero const& zero);
        Matrix(double time, Energies const& eig);
        Matrix(Matrix const&) = delete;
        Matrix(Matrix&&) = delete;
        Matrix& operator=(Matrix const&) = delete;
        Matrix& operator=(Matrix&&) = delete;
        
		int const& I() const { return I_;}
		int const& J() const { return J_;}
        
        ~Matrix();
	protected:
        //int const dim_;        
        int I_, J_;
        data_ptr data_;
		double exponent_;
		
        friend void copyEvolveL(Vector const&, Matrix&, Matrix const&);
		friend void mult(Matrix&, Matrix const&, Matrix const&);
		friend void evolveL(Vector const&, Matrix&);
        
        friend void trace(za::Zahl*, za::Zahl*, Matrix const&);
        friend void accE(double*, za::Zahl, Matrix const&, Energies const&);
        friend void norm(double*, Matrix const&);
        
        friend void axpy(Matrix&, za::Zahl const&, Matrix const&);
        friend void axpy(double*, double, Matrix const&);
	};
	
    void copyEvolveL(Vector const& prop, Matrix& dest, Matrix const& source);
	void mult(Matrix& dest, Matrix const& L, Matrix const& R);
    void evolveL(Vector const& prop, Matrix& arg);
	
    void trace(za::Zahl* Z, za::Zahl* accZ, Matrix const& matrix);
    void accE(double* result, za::Zahl fact, Matrix const& matrix, Energies const& eig);
	void norm(double* norm, Matrix const& matrix);
	
	void axpy(Matrix& dest, za::Zahl const& fact, Matrix const& source);
	void axpy(double* dest, double fact, Matrix const& source);
}


#endif
