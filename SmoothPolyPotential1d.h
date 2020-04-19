/*
 * SmoothPolyPotential1d.h
 *
 *  Created on: Dec 7, 2015
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
#ifndef _SmoothPolyPotential1d_
#define _SmoothPolyPotential1d_

#define _DEFAULT_SOURCE_DIFFERENTIABLITY_ 6

#include <cmath>
#include <vector>
#include "SmoothPolyMollifier1d.h"

class SmoothPolyPotential1d
{

	public:

    SmoothPolyPotential1d()
	{
	initialize();
	}

	SmoothPolyPotential1d(const SmoothPolyPotential1d& S)
	{
    initialize(S);
	}

    SmoothPolyPotential1d(double xPos, double radius, double strength, double laplaceCoeff = 1.0)
	{
    initialize(xPos,radius,strength, laplaceCoeff);
	}

	void initialize()
	{
	xPos              = 0.0;
    radius            = 0.0;
    strength          = 0.0;
    exponent          = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
    laplaceCoeff      = 1.0;

    source.initialize(xPos,radius,strength);
    source.setDifferentiability(this->exponent-1);
	}

    void initialize(const SmoothPolyPotential1d& S)
	{
    strength     = S.strength;
    xPos         = S.xPos;
    radius       = S.radius;
    exponent     = S.exponent;
    laplaceCoeff = S.laplaceCoeff;
    source.initialize(S.source);
	}


	void initialize(double xPos, double radius, double strength, double laplaceCoeff = 1.0)
	{
	this->xPos      = xPos;
    this->radius    = radius;
    this->strength  = strength;
    this->laplaceCoeff      = laplaceCoeff;
    this->exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;

    source.initialize(xPos,radius,strength);
    source.setDifferentiability(this->exponent-1);
	}

	void setRadius(double radius) {this->radius = radius;}
	double getRadius()            {return this->radius;  }

	// 0 <= diffOrder <= 9

	void setSourceDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;

    source.setDifferentiability(this->exponent-1);
	}

	int getSourceDifferentiablity() const
	{
	return this->exponent-1;
	}

    void setLaplaceCoefficient(double laplaceCoeff)
	{
	this->laplaceCoeff = laplaceCoeff;
	}

	double getLaplaceCoefficient() const
	{
	return this->laplaceCoeff;
	}

	SmoothPolyMollifier1d getSource() const
	{
	return source;
	}

	double operator()(double x) const
	{
	double r2       = (x-xPos)*(x-xPos);
	double r2radius = r2/(radius*radius);

	if(r2radius >= 1.0)
	{
	if( x > xPos) return 0.5*(strength/laplaceCoeff)*(x-xPos);
	if( x < xPos) return 0.5*(strength/laplaceCoeff)*(xPos-x);
	}
    return ((strength*radius)/laplaceCoeff)*(evaluation1D_2ndOrder(r2radius));
	}


	//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double)> getEvaluationPtr() const
	{
    std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return std::move(F);
	}
#endif

    double evaluateSource(double x) const
	{
    return source.operator()(x);
	}

	//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double)> getSourceEvaluationPtr() const
	{
    std::function<double(double)> F = [this](double x) {return this->evaluateSource(x);};
	return std::move(F);
	}
#endif


    // Returns the potential due to a unit source with unit radius.
    // r2 is the square of the distance to the center of the source

	double evaluation1D_2ndOrder(double r2)  const
	{
	switch(this->exponent-1)
	{
	case  0 : return    3.0/16.0+(3.0/8.0-r2/16.0)*r2; break;
	case  1 : return    5.0/32.0+(15.0/32.0+(-5.0/32.0+r2/32.0)*r2)*r2; break;
	case  2 : return    35.0/256.0+(35.0/64.0+(-35.0/128.0+(7.0/64.0-5.0/256.0*r2)*r2)*r2)*r2; break;
	case  3 : return    663.0/512.0+(315.0/512.0+(-105.0/256.0+(63.0/256.0+(-45.0/512.0+7.0/512.0*r2)*r2)*r2)*r2)*r2; break;
	case  4 : return    231.0/2048.0+(693.0/1024.0+(-1155.0/2048.0+(231.0/512.0+(-495.0/2048.0+(77.0/1024.0-21.0/2048.0*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  5 : return    429.0/4096.0+(3003.0/4096.0+(-3003.0/4096.0+(3003.0/4096.0+(-2145.0/4096.0+(1001.0/4096.0+(-273.0/4096.0+33.0/4096.0*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  6 : return    6435.0/65536.0+(6435.0/8192.0+(-15015.0/16384.0+(9009.0/8192.0+(-32175.0/32768.0+(5005.0/8192.0+(-4095.0/16384.0+(495.0/8192.0-429.0/65536.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  7 : return    12155.0/131072.0+(109395.0/131072.0+(-36465.0/32768.0+(51051.0/32768.0+(-109395.0/65536.0+(85085.0/65536.0+(-23205.0/32768.0+(8415.0/32768.0+(-7293.0/131072.0+715.0/131072.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
    case  8 : return    46189.0/524288.0+(230945.0/262144.0+(-692835.0/524288.0+(138567.0/65536.0+(-692835.0/262144.0+(323323.0/131072.0+(-440895.0/262144.0+(53295.0/65536.0+(-138567.0/524288.0+(13585.0/262144.0-2431.0/524288.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  9 : return    88179.0/1048576.0+(969969.0/1048576.0+(-1616615.0/1048576.0+(2909907.0/1048576.0+(-2078505.0/524288.0+(2263261.0/524288.0+(-1851759.0/524288.0+(1119195.0/524288.0+(-969969.0/1048576.0+(285285.0/1048576.0+(-51051.0/1048576.0+4199.0/1048576.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	}
    return 0.0;
	}


    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =  xx derivative

	void derivatives(double x, std::vector <double>& derivativeList, int maxOrder = 2) const
   	{
   	switch (maxOrder)
   	{
   	case 0  :  derivativeList.resize(1); break;
   	case 1  :  derivativeList.resize(2); break;
   	case 2  :  derivativeList.resize(3); break;
   	default :  derivativeList.resize(3); break;
   	}

   	double xRadius  = (x-xPos)/radius;
    double r2       = (x-xPos)*(x-xPos);
	double r2radius = r2/(radius*radius);

    int dCount = derivativeList.size();

    if(r2radius > 1.0)
    {
    	if(dCount >= 1)
    	{
    	    if(x > xPos)
    	    {derivativeList[0] =  0.5*(strength/laplaceCoeff)*(x-xPos);}
    		else
    		{derivativeList[0] =  0.5*(strength/laplaceCoeff)*(xPos-x);}
    	}
    	if(dCount >= 2)
    	{
    	    if(x > xPos)
    	    {derivativeList[1] =   0.5*(strength/laplaceCoeff);}
    		else
    		{derivativeList[1] =  -0.5*(strength/laplaceCoeff);}
    	}
    	if(dCount >= 3)
    	{
    		derivativeList[2] =  0.0;
    	}
        return;
    }

    evaluateDerivatives1D(r2radius, xRadius, derivativeList);

    // Scale results

    if(dCount >= 1)
    {
    derivativeList[0] *= ((strength*radius)/laplaceCoeff);
    }

    if(dCount >= 2)
    {
    derivativeList[1] *= (strength/laplaceCoeff);
    }

    if(dCount >= 3)
    {
    derivativeList[2] *= (strength/(radius*laplaceCoeff));
    }

    return;
   	}


	void evaluateDerivatives1D(double r2, double x, std::vector <double>& derivativeList) const
    {
    	int size = derivativeList.size();
    	derivativeList[0] =  evaluation1D_2ndOrder(r2);

    	double dVal = 0.0;
    	if(size >= 2)
    	{

    	// Evaluate derivative of the potential in normalized coordinates, e.g. x = (x-xPos)/radius)

        switch(this->exponent-1)
	    {
	    case  0 : dVal =  (3.0/4.0-r2/4.0)*x; break;
	    case  1 : dVal =  (15.0/16.0+(-5.0/8.0+3.0/16.0*r2)*r2)*x; break;
		case  2 : dVal =  (35.0/32.0+(-35.0/32.0+(21.0/32.0-5.0/32.0*r2)*r2)*r2)*x; break;
		case  3 : dVal =  (315.0/256.0+(-105.0/64.0+(189.0/128.0+(-45.0/64.0+35.0/256.0*r2)*r2)*r2)*r2)*x; break;
		case  4 : dVal =  (693.0/512.0+(-1155.0/512.0+(693.0/256.0+(-495.0/256.0+(385.0/512.0-63.0/512.0*r2)*r2)*r2)*r2)*r2)*x; break;
		case  5 : dVal =  (3003.0/2048.0+(-3003.0/1024.0+(9009.0/2048.0+(-2145.0/512.0+(5005.0/2048.0+(-819.0/1024.0+231.0/2048.0*r2)*r2)*r2)*r2)*r2)*r2)*x; break;
		case  6 : dVal =  (6435.0/4096.0+(-15015.0/4096.0+(27027.0/4096.0+(-32175.0/4096.0+(25025.0/4096.0+(-12285.0/4096.0+(3465.0/4096.0-429.0/4096.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*x; break;
		case  7 : dVal =  (109395.0/65536.0+(-36465.0/8192.0+(153153.0/16384.0+(-109395.0/8192.0+(425425.0/32768.0+(-69615.0/8192.0+(58905.0/16384.0+(-7293.0/8192.0+6435.0/65536.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*x; break;
    	case  8 : dVal =  (230945.0/131072.0+(-692835.0/131072.0+(415701.0/32768.0+(-692835.0/32768.0+(1616615.0/65536.0+(-1322685.0/65536.0+(373065.0/32768.0+(-138567.0/32768.0+(122265.0/131072.0-12155.0/131072.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*x; break;
		case  9 : dVal =  (969969.0/524288.0+(-1616615.0/262144.0+(8729721.0/524288.0+(-2078505.0/65536.0+(11316305.0/262144.0+(-5555277.0/131072.0+(7834365.0/262144.0+(-969969.0/65536.0+(2567565.0/524288.0+(-255255.0/262144.0+46189.0/524288.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*x; break;
		}

		derivativeList[1] = dVal;
		}

    	if(size >= 3)
    	{
    	derivativeList[2] = source.evaluation1D(r2);
    	}
    }



	double     laplaceCoeff;
	double             xPos;
	double           radius;
	double         strength;
	long           exponent;

	SmoothPolyMollifier1d source;

};


#undef _DEFAULT_SOURCE_DIFFERENTIABLITY_

#endif /* _SmoothPolyPotential1d_  */

/* Maple commmands to determine the normalized potential evaluation and derivatives of the normalized potential

restart;  p:=10; g := r -> (969969/524288)*(1-r*r)^p;int(g(r),r=-1..1);
> readlib(C);plot(int(g(r)*(1/2)*std::abs(s-r),r = -1..1),s=-4..4);
> assume(s <= 1);additionally(s >= 0);C(subs(s^2=r2,convert(int(g(r)*(1/2)*std::abs(s-r),r = -1..1),horner)));
> C(subs(s*s=r2,convert(diff(int(g(r)*(1/2)*std::abs(s-r),r = -1..1),s),horner)));

*/

