/*
 * SmoothPolyMollifier1d.h
 *
 *  Created on: June 28,2015
 *      Author: anderson
*/
/*
#############################################################################
#
# Copyright 2015-16-16 Chris Anderson
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Lesser GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# For a copy of the GNU General Public License see
# <http://www.gnu.org/licenses/>.
#
#############################################################################
*/
#include <cmath>
#include <vector>
#include <functional>

#ifndef SMOOTH_POLY_MOLLIFIER_1D_
#define SMOOTH_POLY_MOLLIFIER_1D_

#define DEFAULT_DIFFERENTIABLITY 6

class SmoothPolyMollifier1d
{
	public:

	SmoothPolyMollifier1d()
	{
	initialize();
	}

    SmoothPolyMollifier1d(const SmoothPolyMollifier1d& S)
	{
    initialize(S);
	}

	SmoothPolyMollifier1d(double xPos, double radius, double strength)
	{
    initialize(xPos, radius, strength);
	}

	void initialize()
	{
    radius     = 0.0;
    strength   = 0.0;
    xPos = 0.0;
    exponent     = DEFAULT_DIFFERENTIABLITY + 1 ;
	}

	void initialize(const SmoothPolyMollifier1d& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    radius     = S.radius;
    exponent   = S.exponent;
	}

	void initialize(double xPos, double radius, double strength)
	{
    this->exponent     = DEFAULT_DIFFERENTIABLITY + 1;

    this->strength  = strength;
    this->radius    = radius;
    this->xPos      = xPos;
	}

	void setRadius(double radius) {this->radius = radius;}
	double getRadius()            {return this->radius;  }

    // 0 <= diffOrder <= 9

	void setDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;
	}

	int getDifferentiability() const
	{
	return this->exponent-1;
	}

	double operator()(double x) const
	{
	double r2radius = ((x-xPos)*(x-xPos))/(radius*radius);
    if(r2radius >= 1.0) return 0.0;
    return (strength/radius)*evaluation1D(r2radius);
	}

    //
    // Normalized evaluation routine:
    //
    // Input:
    //
    // A value r2 = r*r for r in [0,1]
    //
    // Output:
    //
    // Returns sigma*(1-r2)^(exponent)
    //
    // where sigma is the pre-factor so that the  mollifier has unit integral over [-1,1].
    //

    double evaluation1D(double r2) const
    {
    double u = 1.0-r2;
    switch(this->exponent-1)
    {
	case  0 : return(3.0/4.0)*u; break;
    case  1 : return(15.0/16.0)*u*u; break;
    case  2 : return(35.0/32.0)*u*u*u; break;
    case  3 : return(315.0/256.0)*u*u*u*u; break;
    case  4 : return(693.0/512.0)*u*u*u*u*u; break;
    case  5 : return(3003.0/2048.0)*u*u*u*u*u*u; break;
    case  6 : return(6435.0/4096.0)*u*u*u*u*u*u*u; break;
    case  7 : return(109395.0/65536.0)*u*u*u*u*u*u*u*u; break;
    case  8 : return(230945.0/131072.0)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return(969969.0/524288.0)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
    }

    //  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return F;
	}

//  Returns a std::function that is bound to the derivative operator of *this


	std::function<double(double)> getDerivativeEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->derivative(x);};
	return F;
	}


   	double derivative(double x) const
   	{
   	double xRadius  = (x-xPos)/radius;
    double r2radius = xRadius*xRadius;

    double derivativeList[2];

    if(r2radius >= 1.0) {return 0.0;}

    evaluateDerivatives1D(r2radius, xRadius, &derivativeList[0], 2);
    return derivativeList[1]*(strength/(radius*radius));
    }

    // Derivative ordering:
    // derivativeList[0] = 0th  derivative
    // derivativeList[1] =   x  derivative
    // derivativeList[2] =  xx  derivative
    // derivativeList[3] = xxx derivative

   	void derivatives(double x, std::vector <double>& derivativeList, int maxOrder = 2) const
   	{
   	switch (maxOrder)
   	{
   	case 0  :  derivativeList.resize(1); break;
   	case 1  :  derivativeList.resize(2); break;
   	case 2  :  derivativeList.resize(3); break;
   	case 3  :  derivativeList.resize(4); break;
   	default :  derivativeList.resize(3); break;
   	}

   	derivativeList.resize(maxOrder+1);

   	double xRadius  = (x-xPos)/radius;
    double r2radius = xRadius*xRadius;

    int dCount = (int)derivativeList.size();

    if(r2radius >= 1.0)
    {
    for(long i = 0; i < dCount; i++)
    {
    derivativeList[i] = 0.0;
    }
    return;
    }

    evaluateDerivatives1D(r2radius, xRadius, &derivativeList[0],dCount);

    if(dCount >= 1)
    {
    derivativeList[0] *= (strength/radius);
    }

    if(dCount >= 2)
    {
	derivativeList[1] *= (strength/(radius*radius));
    }

    if(dCount >= 3)
    {
	derivativeList[2] *= (strength/(radius*radius*radius));
    }

    if(dCount >= 4)
    {
	derivativeList[3] *= (strength/(radius*radius*radius*radius));
    }

    return;
   	}

   	private :

	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	//    Utility functions for evaluating the derivatives of the mollifer
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	// The derivative for r2 <= 1  in 1D

    // Derivative ordering:
    // derivativeList[0] = 0th  derivative
    // derivativeList[1] =   x  derivative
    // derivativeList[2] =  xx  derivative
    // derivativeList[3] =  xxz derivative


    // 0 <=  differentiability order <= 9

    // u            = sigma*(1-r^2)^M
	// du/dx        = sigma*M*(1-r^2)^(M-1)*(-2 r) (x/r) = -2.0*sigma*M*(1-r^2)^(M-1)*x
	// d^2u/dx^2    = 4.0*sigma*M*(M-1)*(1-r^2)^(M-2)*x*x - 2.0*sigma*M*(1-r^2)^(M-1)
	// d^3u/dx^3    =-8.0*sigma*M*(M-1)*(M-2)*uMM3*x*x*x + 12.0*sigma*M*(M-1)*uMM2*x;

    void evaluateDerivatives1D(double r2, double x, double* derivativeList, int size) const
    {
    	if(size == 1)  {derivativeList[0] = evaluation1D(r2); return;}

    	double sigma = 0.0;

		switch(this->exponent-1)
		{
		case  0 : sigma = (3.0/4.0); break;
		case  1 : sigma = (15.0/16.0); break;
		case  2 : sigma = (35.0/32.0); break;
		case  3 : sigma = (315.0/256.0); break;
		case  4 : sigma = (693.0/512.0); break;
		case  5 : sigma = (3003.0/2048.0); break;
		case  6 : sigma = (6435.0/4096.0); break;
    	case  7 : sigma = (109395.0/65536.0); break;
    	case  8 : sigma = (230945.0/131072.0); break;
    	case  9 : sigma = (969969.0/524288.0); break;
		}

		double u     = 1.0-r2;

    	double uMM1 = 0.0; // (1-r^2)^(exponent-1)
    	double uMM2 = 0.0; // (1-r^2)^(exponent-2)
    	double uMM3 = 0.0; // (1-r^2)^(exponent-3)
    	long M = exponent;


    	switch(this->exponent-3)
    	{
		case  0 : uMM3 = 1.0; break;
    	case  1 : uMM3 = u; break;
    	case  2 : uMM3 = u*u; break;
    	case  3 : uMM3 = u*u*u; break;
    	case  4 : uMM3 = u*u*u*u; break;
    	case  5 : uMM3 = u*u*u*u*u; break;
    	case  6 : uMM3 = u*u*u*u*u*u; break;
    	case  7 : uMM3 = u*u*u*u*u*u*u; break;
    	default : uMM3 = 0.0; break;
    	}

    	uMM2 = uMM3*u;
    	uMM1 = uMM2*u;

    	if(exponent-3 == -2) {uMM3 = 0.0; uMM2 = 0.0; uMM1 = 1.0;}
    	if(exponent-3 == -1) {uMM3 = 0.0; uMM2 = 1.0; uMM1 =   u;}

    	if(size >= 2)
    	{
    	derivativeList[0] =  sigma*uMM1*u;
    	derivativeList[1] = -2.0*sigma*uMM1*M*x;
    	}

    	if(size >= 3)
    	{
    	derivativeList[2] =  4.0*sigma*M*(M-1)*uMM2*x*x - 2.0*sigma*uMM1*M;
    	}

    	if(size >= 4)
    	{
        derivativeList[3] = -8.0*sigma*M*(M-1)*(M-2)*uMM3*x*x*x + 12.0*sigma*M*(M-1)*uMM2*x;
        }

}

    double    strength;    // The strength of the mollifier
    double        xPos;    // The position of the mollifier

    double        radius;   // The radius of the mollifier
    int         exponent;   // The exponent of the mollifier
};


#undef DEFAULT_DIFFERENTIABLITY
#endif /* SMOOTHPOLYMOLLIFIER_H_ */
