/*
 * SmoothPolyMollifier2d.h
 *
 *  Created on: Feb 13, 2014
 *      Author: anderson
*/
/*
#############################################################################
#
# Copyright 2014-2015 Chris Anderson
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

#ifndef _SmoothPolyMollifier2d_
#define _SmoothPolyMollifier2d_

#define _DEFAULT_DIFFERENTIABLITY_ 6

class SmoothPolyMollifier2d
{
	public:

	SmoothPolyMollifier2d()
	{
	initialize();
	}

    SmoothPolyMollifier2d(const SmoothPolyMollifier2d& S)
	{
    initialize(S);
	}

    SmoothPolyMollifier2d(double xPos, double yPos, double radius, double strength)
	{
    initialize(xPos,yPos,radius,strength);
	}

	void initialize()
	{
    radius     = 0.0;

    strength  = 0.0;
    xPos = 0.0; yPos = 0.0;
    exponent     = _DEFAULT_DIFFERENTIABLITY_ + 1 ;
	}

	void initialize(const SmoothPolyMollifier2d& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    yPos       = S.yPos;
    radius     = S.radius;
    exponent   = S.exponent;
	}

	void initialize(double xPos, double yPos, double radius, double strength)
	{
    this->exponent     = _DEFAULT_DIFFERENTIABLITY_ + 1;

    this->strength  = strength;
    this->radius    = radius;

    this->xPos      = xPos;
    this->yPos      = yPos;
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


	double operator()(double x, double y) const
	{
    double r2radius = ((x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos))/(radius*radius);
    if(r2radius >= 1.0) return 0.0;
    return (strength/(radius*radius))*evaluation2D(r2radius);
	}

    //
    // Normalized evaluation routine:
    //
    // Input:  A value r2 = r*r for r in [0,1]
    //
    // Output:  Returns (sigma/(2*Pi))*(1-r2)^(exponent)
    // where sigma is the prefactor so that the mollifier
    // has unit integral over the unit ball.
    //

	double evaluation2D(double r2) const
	{
    double u = 1.0-r2;
    double oneOver2pi = 0.159154943091895336;

	switch(this->exponent-1)
    {
	case  0 : return   oneOver2pi*2*(exponent+1)*u;  break;
    case  1 : return   oneOver2pi*2*(exponent+1)*u*u; break;
    case  2 : return   oneOver2pi*2*(exponent+1)*u*u*u; break;
    case  3 : return   oneOver2pi*2*(exponent+1)*u*u*u*u; break;
    case  4 : return   oneOver2pi*2*(exponent+1)*u*u*u*u*u; break;
    case  5 : return   oneOver2pi*2*(exponent+1)*u*u*u*u*u*u; break;
    case  6 : return   oneOver2pi*2*(exponent+1)*u*u*u*u*u*u*u; break;
    case  7 : return   oneOver2pi*2*(exponent+1)*u*u*u*u*u*u*u*u; break;
    case  8 : return   oneOver2pi*2*(exponent+1)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return   oneOver2pi*2*(exponent+1)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
	}

	//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double,double)> getEvaluationPtr() const
	{
    std::function<double(double,double)> F = [this](double x,double y) {return this->operator()(x,y);};
	return std::move(F);
	}


    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =  xx derivative
    // derivativeList[4] =  xy derivative
    // derivativeList[5] =  yy derivative

   	void derivatives(double x, double y, std::vector <double>& derivativeList,int maxOrder = 2) const
   	{
   	switch (maxOrder)
   	{
   	case 0  :  derivativeList.resize(1); break;
   	case 1  :  derivativeList.resize(3); break;
   	case 2  :  derivativeList.resize(6); break;
   	default :  derivativeList.resize(6); break;
   	}


   	double xRadius    = (x-xPos)/radius;
   	double yRadius    = (y-yPos)/radius;
    double r2radius   = xRadius*xRadius + yRadius*yRadius;

    double oneOver2pi = 0.159154943091895336;

    int dCount = derivativeList.size();

    if(r2radius >= 1.0)
    {
    for(long i = 0; i < dCount; i++)
    {
    derivativeList[i] = 0.0;
    }
    return;
    }

    evaluateDerivatives2D(r2radius, xRadius, yRadius, derivativeList);

    if(dCount >= 1)
    {
    derivativeList[0] *= oneOver2pi*(strength/(radius*radius));
    }

    if(dCount >= 3)
    {
    	for(long i = 1; i <= 2; i++)
    	{
    	derivativeList[i] *= oneOver2pi*(strength/(radius*radius*radius));
    	}
    }

    if(dCount >= 6)
    {
    	for(long i = 3; i <= 5; i++)
    	{
    	derivativeList[i] *= oneOver2pi*(strength/(radius*radius*radius*radius));
    	}
    }
    return;
   	}

   	// Internal normalized routines

   	private:





	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	//    Utility functions for evaluating the derivatives of the mollifer
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Creates 2*Pi* value of the derivative of the normalized mollifier.


    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =  xx derivative
    // derivativeList[4] =  xy derivative
    // derivativeList[5] =  yy derivative
    //
    // 0 <=  differentiability order <= 9
    //
	// u(r)         = sigma*(1-r^2)^M
	// du/dx(r)     = sigma*M*(1-r^2)^(M-1)*(-2 r) (x/r) = -2.0*sigma*M*(1-r^2)^(M-1)*x
	// d^2u/dx^2(r) = 4.0*sigma*M*(M-1)*(1-r^2)^(M-2)*x*x - 2.0*sigma*M*(1-r^2)^(M-1)
    //
    void evaluateDerivatives2D(double r2, double x, double y, std::vector <double>& derivativeList) const
    {
    	int size = derivativeList.size();

    	if(size == 1)
    	{derivativeList[0] = evaluation2D(r2)*6.28318530717958648; return;}

    	double sigma = 2.0*(exponent+1);
		double u     = 1.0-r2;

    	double uMM1 = 1.0; // (1-r^2)^(exponent-1)
    	double uMM2 = 1.0; // (1-r^2)^(exponent-2)
    	double M = exponent;


    	switch(this->exponent-2)
    	{
		case  0 : uMM2 = 1.0;  break;
    	case  1 : uMM2 = u; break;
    	case  2 : uMM2 = u*u; break;
    	case  3 : uMM2 = u*u*u; break;
    	case  4 : uMM2 = u*u*u*u; break;
    	case  5 : uMM2 = u*u*u*u*u; break;
    	case  6 : uMM2 = u*u*u*u*u*u; break;
    	case  7 : uMM2 = u*u*u*u*u*u*u; break;
    	case  8 : uMM2 = u*u*u*u*u*u*u*u; break;
    	}

    	uMM1 = uMM2*u;
        uMM1 = uMM2*u;
    	if(exponent-2 == -1)
    	{
    	uMM2 = 0.0;
    	uMM1 = 1.0;
    	}

    	if(size >= 3)
    	{
    	derivativeList[0] =  sigma*uMM1*u;
    	derivativeList[1] = -2.0*sigma*uMM1*M*x;
    	derivativeList[2] = -2.0*sigma*uMM1*M*y;
    	}

    	double MM2factor = 4.0*sigma*M*(M-1.0);
    	if(size >= 6)
    	{
    	derivativeList[3] = MM2factor*uMM2*x*x - 2.0*sigma*uMM1*M;
    	derivativeList[4] = MM2factor*uMM2*x*y;
    	derivativeList[5] = MM2factor*uMM2*y*y - 2.0*sigma*uMM1*M;
    	}
    }



    double    strength;    // The strength of the mollifier
    double        xPos;    // The position of the mollifier
    double        yPos;

    double        radius;   // The radius of the mollifier
    int         exponent;   // The exponent of the mollifier
};


#undef _DEFAULT_DIFFERENTIABLITY_
#endif /* SMOOTHPOLYMOLLIFIER_H_ */
