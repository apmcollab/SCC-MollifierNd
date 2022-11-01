/*
 * SmoothPolyStep.h
 *
 *  Created on: Feb 13, 2014
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
#include <functional>

//
// Class SmoothPolyStep
//
// The operator()(x) member functions returns the value of a smooth step
// function that increases monotonically from a value 0 to a value strength
// as x goes from (xPos - delta/2) to (xPos + delta/2) where
// delta is the width of the step.
//
//
//                                                 | strength
//
//    0 |<--    transition from 0 to strength  --> |
//
//      |--------------------|---------------------|
// xPos - delta/2              xPos           xPos + delta/2
//
//
// The null initializer sets xPos = 0.0 and strength = 1 so that it is
// a mollified Heaviside function.
//
// This function is symmetric about x = xPos.
//
// If the stepDownFlag is set to true, then the function transitions
// from strength down to 0 instead of from 0 to strength.
//
// The differentiability of the step function can set, and ranges
// from 1 to 10 with the default being 6
//
// These functions are created by scaling and translating the
// integral of the smooth polynomial bump function B(x)
// from x = -1 to s for s <= 1 where
//
// B(x) = sigma*(1-x^2)^i       for x in [-1,1]
//      =     0                 elsewhere
//
// for integer i in [1,10].
//
// The normalization constants for B(x) and a description of the origin of he mollifiers
// is given in
//
// C.R. Anderson,”Compact Polynomial Mollifiers For Poisson's Equation”,
// Technical Report CAM-14-43, UCLA.
//
// Version 6/28/2015
//
#include "SmoothPolyMollifier1d.h"

#ifndef SMOOTH_POLY_STEP_
#define SMOOTH_POLY_STEP_

#define DEFAULT_DIFFERENTIABLITY 6

class SmoothPolyStep 
{
	public:

	SmoothPolyStep()
	{
	initialize();
	}

	SmoothPolyStep(const SmoothPolyStep& S)
	{
    strength     = S.strength;
    xPos         = S.xPos;
    deltaStar    = S.deltaStar;
    exponent     = S.exponent;
    stepDownFlag = S.stepDownFlag;
	}

	SmoothPolyStep(double xPos,double delta, double strength, bool stepDownFlag = false)
	{
    initialize(xPos, delta, strength, stepDownFlag);
	}

	void initialize()
	{
    deltaStar    = 0.0;
    exponent     = DEFAULT_DIFFERENTIABLITY;
    strength     = 1.0;
    xPos         = 0.0;
    stepDownFlag = false;
	}

	void initialize(const SmoothPolyStep& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    deltaStar  = S.deltaStar;
    exponent   = S.exponent;
    stepDownFlag = S.stepDownFlag;
	}


	void initialize(double xPos,double delta, double strength, bool stepDownFlag = false)
	{
    this->exponent     = DEFAULT_DIFFERENTIABLITY;
    this->strength     = strength;
    this->deltaStar    = delta/2.0;
    this->xPos         = xPos;
    this->stepDownFlag = stepDownFlag;
	}

	void setTransitionDistance(double delta) {this->deltaStar  = delta/2.0;}
	double getTransitionDistance()           {return   this->deltaStar*2.0;}

	// 1 <= diffOrder <= 10

	void setDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder;
	if(exponent > 10) {this->exponent = 10;}
    if(exponent < 1 ) {this->exponent = 1;}
	}

	int getDifferentiablity() const
	{
	return this->exponent;
	}

	void setStepDownFlag(bool stepDownFlag = true)
	{
	this->stepDownFlag =  stepDownFlag;
	}

	void clearStepDownFlag()
	{
	this->stepDownFlag =  false;
	}

	// Returns value of transition function
	// centered at xPos and of width delta/2.0


	double operator()(double x) const
	{
    double xDelta = (x-xPos)/deltaStar;

    // Translate, scale and then call normalized step function

    if(stepDownFlag)
    {
    if(xDelta <= -1.0) {return strength;}
    if(xDelta >=  1.0) {return 0.0;}
    return strength*(1.0-evaluate(xDelta));
    }

    if(xDelta <= -1.0) {return 0.0;}
    if(xDelta >=  1.0) {return strength;}
    return strength*evaluate(xDelta);
	}

	// Returns value of the transition function
	// centered at xPos and of width delta/2.0
	// The value of delta associated with this
	// instance is NOT set by the parameter delta

	double  operator()(double x, double delta) const
	{

	// Translate, scale and then call normalized step function

	double deltaStar = delta/2.0;
	double xDelta = (x-xPos)/deltaStar;

    if(stepDownFlag)
    {
    if(xDelta <= -1.0) {return strength;}
    if(xDelta >=  1.0) {return 0.0;}
    return strength*(1.0-evaluate(xDelta));
    }

	if(xDelta <= -1.0) {return 0.0;}
    if(xDelta >=  1.0) {return strength;}
	return strength*evaluate(xDelta);
	}

	// Normalized step function
	//
	// Returns the value of the integral from -1 to s of alpha*(1 - s^2)^p where
	// p is the exponent (= differentiability)
	// and alpha is a normalization factor so that
	// the value 0 is returned at s = 0.0 and
	// the value 1 is returned at s = 1.0
	//
	//
    double evaluate(double s) const
    {
    double t1 = s*s;
	double rVal = 0.0;
    switch(this->exponent-1)
    {
    case  0 :rVal = 1.0/2.0+(3.0/4.0-t1/4)*s; break;
    case  1 :rVal = 1.0/2.0+(15.0/16.0+(-5.0/8.0+3.0/16.0*t1)*t1)*s;  break;
    case  2 :rVal = 1.0/2.0+(35.0/32.0+(-35.0/32.0+(21.0/32.0-5.0/32.0*t1)*t1)*t1)*s;  break;
    case  3 :rVal = 1.0/2.0+(315.0/256.0+(-105.0/64.0+(189.0/128.0+(-45.0/64.0+35.0/256.0*t1)*t1)*t1)*t1)*s;  break;
    case  4 :rVal = 1.0/2.0+(693.0/512.0+(-1155.0/512.0+(693.0/256.0+(-495.0/256.0+(385.0/512.0-63.0/512.0*t1)*t1)*t1)*t1)*t1)*s;  break;
    case  5 :rVal = 1.0/2.0+(3003.0/2048.0+(-3003.0/1024.0+(9009.0/2048.0+(-2145.0/512.0+(5005.0/2048.0+(-819.0/1024.0+231.0/2048.0*t1)*t1)*t1)*t1)*t1)*t1)*s; break;
    case  6 :rVal = 1.0/2.0+(6435.0/4096.0+(-15015.0/4096.0+(27027.0/4096.0+(-32175.0/4096.0+(25025.0/4096.0+(-12285.0/4096.0+(3465.0/4096.0-429.0/4096.0*t1)*t1)*t1)*t1)*t1)*t1)*t1)*s; break;
    case  7 :rVal = 1.0/2.0+(109395.0/65536.0+(-36465.0/8192.0+(153153.0/16384.0+(-109395.0/8192.0+(425425.0/32768.0+(-69615.0/8192.0+(58905.0/16384.0+(-7293.0/8192.0+6435.0/65536.0*t1)*t1)*t1)*t1)*t1)*t1)*t1)*t1)*s; break;
    case  8 :rVal = 1.0/2.0+(230945.0/131072.0+(-692835.0/131072.0+(415701.0/32768.0+(-692835.0/32768.0+(1616615.0/65536.0+(-1322685.0/65536.0+(373065.0/32768.0+(-138567.0/32768.0+(122265.0/131072.0-12155.0/131072.0*t1)*t1)*t1)*t1)*t1)*t1)*t1)*t1)*t1)*s; break;
    case  9 :rVal = 1.0/2.0+(969969.0/524288.0+(-1616615.0/262144.0+(8729721.0/524288.0+(-2078505.0/65536.0+(11316305.0/262144.0+(-5555277.0/131072.0+(7834365.0/262144.0+(-969969.0/65536.0+(2567565.0/524288.0+(-255255.0/262144.0+46189.0/524288.0*t1)*t1)*t1)*t1)*t1)*t1)*t1)*t1)*t1)*t1)*s; break;
    }
    return rVal;
    }

    //  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double)> getEvaluationPtr()
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return F;
	}


    // Returns the smooth polynomial mollifier upon which this step function is based

	SmoothPolyMollifier1d getDerivative() const
	{
	SmoothPolyMollifier1d M;

	if(stepDownFlag)
	{M.initialize(xPos,deltaStar,-strength);}
	else
	{M.initialize(xPos,deltaStar,strength);}

    M.setDifferentiability(this->exponent-1);
    return M;
	}

    double    strength;    // The strength of the step
    double        xPos;    // The position of the step
    double    deltaStar;   // The width of the step
    int         exponent;  // The exponent of the mollifier function upon which the smooth step function is based 1 <= exponent <= 10;
    bool    stepDownFlag;  // Evaluation returns transition from strength to 0 instead of 0 to strength
};


#undef DEFAULT_DIFFERENTIABLITY

#endif /* SMOOTHPOLYMOLLIFIER_H_ */
