/*
 * HighOrderPolyPotential2d.h
 *
 *  Created on: Feb 16, 2014
 *      Author: anderson
 *
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
using namespace std;

#ifndef _HighOrderPolyPotential2d_
#define _HighOrderPolyPotential2d_

#define _DEFAULT_ORDER_      4
#define _DEFAULT_SOURCE_DIFFERENTIABLITY_ 6

class HighOrderPolyPotential2d
{
	public:

	HighOrderPolyPotential2d()
	{
	initialize();
	}

	HighOrderPolyPotential2d(const HighOrderPolyPotential2d& S)
	{
    xPos       = S.xPos;
    yPos       = S.yPos;
    radius       = S.radius;
    laplaceCoeff = S.laplaceCoeff;
    strength     = S.strength;

    order      = S.order;
    exponent   = S.exponent;
    logRadius  = S.logRadius;
	}

    HighOrderPolyPotential2d(double xPos, double yPos, double radius,double strength, double laplaceCoeff = 1.0)
	{
    initialize(xPos,yPos,radius,strength,laplaceCoeff);
	}

	void initialize()
	{
	xPos      = 0.0;
    yPos      = 0.0;
    radius         = 0.0;
    laplaceCoeff   = 1.0;
    strength       = 0.0;

    order     = _DEFAULT_ORDER_;
    exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
    logRadius = 0.0;
	}

	void initialize(double xPos, double yPos, double radius, double strength, double laplaceCoeff = 1.0)
	{
    this->xPos      = xPos;
    this->yPos      = yPos;
    this->radius    = radius;
    this->laplaceCoeff  = laplaceCoeff;
    this->strength  = strength;

    this->exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
    this->order     = _DEFAULT_ORDER_;
    this->logRadius = log(radius);
	}

	void setRadius(double radius) {this->radius = radius;}
	double getRadius()            {return this->radius;  }

	void setOrder(int order)
	{
	this->order = order;
	if((this->order != 2)&&(this->order != 4)&&(this->order != 6))
	{
	this->order = _DEFAULT_ORDER_;
	}
	}

	int getOrder() const
	{return this->order;}



	// 0 <= diffOrder <= 9

	void setSourceDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;
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

	// In 2D a shift is added since the mollified version of the potential is phi(r/radius)
	// where phi(r) = u(r/radius) r/radius <= 1 and  1/(2pi) log (r/radius) for r/radius> 1.
	//
	// To match the asymptotic behavior 1/(2pi) log(r) necessitates adding the
	// constant value 1/(2*pi)log(radius).


	double operator()(double x, double y) const
	{
	double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos);
	double r2radius = r2/(radius*radius);
	if(r2radius >= 1.0) {return (strength/laplaceCoeff)*.15915494309190*log(sqrt(r2));}
    switch(order)
    {
    case 2 : return ((strength/laplaceCoeff)*.15915494309190)*(evaluation2D_2ndOrder(r2radius) + logRadius); break;
    case 4 : return ((strength/laplaceCoeff)*.15915494309190)*(evaluation2D_4thOrder(r2radius) + logRadius); break;
    case 6 : return ((strength/laplaceCoeff)*.15915494309190)*(evaluation2D_6thOrder(r2radius) + logRadius); break;
    }
    return 0.0;
	}

//  Returns a std::function that is bound to the evaluation operator of *this

#if __cplusplus > 199711L
	std::function<double(double,double)> getEvaluationPtr() const
	{
	std::function<double(double,double)> F = [this](double x,double y) {return this->operator()(x,y);};
	return std::move(F);
	}
#endif

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative

    // 0 <=  differentiability order <= 9

   	void derivatives(double x, double y, vector <double>& derivativeList, int maxOrder = 1) const
   	{
   	switch (maxOrder)
   	{
   	case 0   :  derivativeList.resize(1); break;
   	case 1   :  derivativeList.resize(3); break;
   	default  :  derivativeList.resize(3); break;
   	}

   	double xRadius  = (x-xPos)/radius;
   	double yRadius  = (y-yPos)/radius;
    double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos);
	double r2radius = r2/(radius*radius);
	double r;

    int dCount = derivativeList.size();

    if(r2radius <= 1.0)
    {
    evaluateDerivatives2D(r2radius, xRadius, yRadius, derivativeList);
    }
    else
    {
    	r = sqrt(r2);
    	if(dCount >= 1)
    	{
    		derivativeList[0] =  (strength/laplaceCoeff)*.15915494309190*log(r);
    	}
    	if(dCount >= 3)
    	{
    		derivativeList[1] =  (strength/laplaceCoeff)*.15915494309190*((x-xPos)/r2);
    		derivativeList[2] =  (strength/laplaceCoeff)*.15915494309190*((y-yPos)/r2);
    	}
    return;
    }

    // Scale results

    if(dCount >= 1)
    {
    derivativeList[0] += logRadius;
    derivativeList[0] *= ((strength/laplaceCoeff)*.15915494309190);
    }

    if(dCount >= 3)
    {
    	for(long i = 1; i <= 2; i++)
    	{
    	derivativeList[i] *= (strength/(laplaceCoeff*radius*6.28318530717958648));
    	}
    }

    return;
   	}


//  Internal normalized evaluation routines

private :

    // returns 2*pi*mollifier value

	double evaluation2D_2ndOrder(double r2) const
	{
	switch(this->exponent-1)
	{
	case  0 : return    2*(exponent+1)*((1.0/4.0-r2/16)*r2 - 3.0/16.0); break;
	case  1 : return    2*(exponent+1)*((1.0/4.0+(-1.0/8.0+r2/36.0)*r2)*r2 - 11.0/72.0); break;
	case  2 : return    2*(exponent+1)*((1.0/4.0+(-3.0/16.0+(1.0/12.0-r2/64.0)*r2)*r2)*r2 -  25.0/192.0); break;
	case  3 : return    2*(exponent+1)*((1.0/4.0+(-1.0/4.0+(1.0/6.0+(-1.0/16.0+r2/100.0)*r2)*r2)*r2)*r2 - 137.0/1200.0); break;
	case  4 : return    2*(exponent+1)*((1.0/4.0+(-5.0/16.0+(5.0/18.0+(-5.0/32.0+(1.0/20.0-r2/144.0)*r2)*r2)*r2)*r2)*r2 - 49.0/480.0); break;
	case  5 : return    2*(exponent+1)*((1.0/4.0+(-3.0/8.0+(5.0/12.0+(-5.0/16.0+(3.0/20.0+(-1.0/24.0+r2/196)*r2)*r2)*r2)*r2)*r2)*r2 - 363.0/3920.0); break;
	case  6 : return    2*(exponent+1)*((1.0/4.0+(-7.0/16.0+(7.0/12.0+(-35.0/64.0+(7.0/20.0+(-7.0/48.0+(1.0/28.0-r2/256.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2 -  761.0/8960.0); break;
	case  7 : return    2*(exponent+1)*((1.0/4.0+(-1.0/2.0+(7.0/9.0+(-7.0/8.0+(7.0/10.0+(-7.0/18.0+(1.0/7.0+(-1.0/32.0+r2/324.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 7129.0/90720.0); break;
	case  8 : return    2*(exponent+1)*((1.0/4.0+(-9.0/16.0+(1.0+(-21.0/16.0+(63.0/50.0+(-7.0/8.0+(3.0/7.0+(-9.0/64.0+(1.0/36.0-r2/400.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 7381.0/100800.0); break;
	case  9 : return    2*(exponent+1)*((1.0/4.0+(-5.0/8.0+(5.0/4.0+(-15.0/8.0+(21.0/10.0+(-7.0/4.0+(15.0/14.0+(-15.0/32.0+(5.0/36.0+(-1.0/40.0+r2/484.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 83711.0/1219680.0);break;
	}
    return 0.0;
	}

    // returns 2*pi* value

	double evaluation2D_4thOrder(double r2) const
	{
	switch(this->exponent-1)
    {
	case  0 : return (-12.0)*(17.0/144.0+(-1.0/4.0+(3.0/16.0-r2/18)*r2)*r2);           break;
    case  1 : return (-24.0)*(37.0/576.0+(-1.0/6.0+(3.0/16.0+(-1.0/9.0+5.0/192.0*r2)*r2)*r2)*r2);   break;
    case  2 : return (-40.0)*(197.0/4800.0+(-1.0/8.0+(3.0/16.0+(-1.0/6.0+(5.0/64.0-3.0/200.0*r2)*r2)*r2)*r2)*r2);  break;
    case  3 : return(-60.0)*(23.0/800.0+(-1.0/10.0+(3.0/16.0+(-2.0/9.0+(5.0/32.0+(-3.0/50.0+7.0/720.0*r2)*r2)*r2)*r2)*r2)*r2);   break;
    case  4 : return(-84.0)*( 503.0/23520.0+(-1.0/12.0+(3.0/16.0+(-5.0/18.0+(25.0/96.0+(-3.0/20.0+(7.0/144.0-r2/147)*r2)*r2)*r2)*r2)*r2)*r2);    break;
    case  5 : return(-112.0)*(1041.0/62720.0+(-1.0/14.0+(3.0/16.0+(-1.0/3.0+(25.0/64.0+(-3.0/10.0+(7.0/48.0+(-2.0/49.0+9.0/1792.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);   break;
    case  6 : return(-144.0)*(9649.0/725760.0+(-1.0/16.0+(3.0/16.0+(-7.0/18.0+(35.0/64.0+(-21.0/40.0+(49.0/144.0+(-1.0/7.0+(9.0/256.0-5.0/1296.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  7 : return(-180.0)*(9901.0/907200.0+(-1.0/18.0+(3.0/16.0+(-4.0/9.0+(35.0/48.0+(-21.0/25.0+(49.0/72.0+(-8.0/21.0+(9.0/64.0+(-5.0/162.0+11.0/3600.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  8 : return(-220.0)*(111431.0/12196800.0+(-1.0/20.0+(3.0/16.0+(-1.0/2.0+(15.0/16.0+(-63.0/50.0+(49.0/40.0+(-6.0/7.0+(27.0/64.0+(-5.0/36.0+(11.0/400.0-3.0/1210.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);    break;
    case  9 : return(-264.0)*(113741.0/14636160.0+(-1.0/22.0+(3.0/16.0+(-5.0/9.0+(75.0/64.0+(-9.0/5.0+(49.0/24.0+(-12.0/7.0+(135.0/128.0+(-25.0/54.0+(11.0/80.0+(-3.0/121.0+13.0/6336.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    }
    return 0.0;
	}

    // returns 2*pi*value

    double evaluation2D_6thOrder(double r2) const
	{
   	switch(this->exponent-1)
    {
	case  0 : return(24.0)*(-43.0/576.0+(1.0/4.0+(-3.0/8.0+(5.0/18.0-5.0/64.0*r2)*r2)*r2)*r2);  break;
    case  1 : return(60.0)*(-227.0/7200.0+(1.0/8.0+(-1.0/4.0+(5.0/18.0+(-5.0/32.0+7.0/200.0*r2)*r2)*r2)*r2)*r2);  break;
    case  2 : return(120.0)*(-79.0/4800.0+(3.0/40.0+(-3.0/16.0+(5.0/18.0+(-15.0/64.0+(21.0/200.0-7.0/360.0*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  3 : return(210.0)*(-191.0/19600.0+(1.0/20.0+(-3.0/20.0+(5.0/18.0+(-5.0/16.0+(21.0/100.0+(-7.0/90.0+3.0/245.0*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  4 : return(336.0)*(-1181.0/188160.0+(1.0/28.0+(-1.0/8.0+(5.0/18.0+(-25.0/64.0+(7.0/20.0+(-7.0/36.0+(3.0/49.0-15.0/1792.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  5 : return(504.0)*(-10909.0/2540160.0+(3.0/112.0+(-3.0/28.0+(5.0/18.0+(-15.0/32.0+(21.0/40.0+(-7.0/18.0+(9.0/49.0+(-45.0/896.0+55.0/9072.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  6 : return(720.0)*(-11161.0/3628800.0+(1.0/48.0+(-3.0/32.0+(5.0/18.0+(-35.0/64.0+(147.0/200.0+(-49.0/72.0+(3.0/7.0+(-45.0/256.0+(55.0/1296.0-11.0/2400.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
    case  7 : return(990.0)*(-125291.0/54885600.0+(1.0/60.0+(-1.0/12.0+(5.0/18.0+(-5.0/8.0+(49.0/50.0+(-49.0/45.0+(6.0/7.0+(-15.0/32.0+(55.0/324.0+(-11.0/300.0+13.0/3630.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  8 : return(1320.0)*(-127601.0/73180800.0+(3.0/220.0+(-3.0/40.0+(5.0/18.0+(-45.0/64.0+(63.0/50.0+(-49.0/30.0+(54.0/35.0+(-135.0/128.0+(55.0/108.0+(-33.0/200.0+(39.0/1210.0-91.0/31680.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  9 : return(1716.0)*(-1686533.0/1236755520.0+(1.0/88.0+(-3.0/44.0+(5.0/18.0+(-25.0/32.0+(63.0/40.0+(-7.0/3.0+(18.0/7.0+(-135.0/64.0+(275.0/216.0+(-11.0/20.0+(39.0/242.0+(-91.0/3168.0+35.0/14872.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    }
    return 0.0;
	}

	//
	// Routines required for the evaluation of the derivatives of the potentials
	//

    // Returns 2*Pi*2*q'(r^2) where u(r) = q(r^2)

	double evaluate2xDq2D_2ndOrder(double r2) const
	{

	// Switch based on the differentiability of the underlying mollifer

	switch(this->exponent-1)
	{
	case  0 : return     2*(exponent+1)*(1.0/2.0-r2/4); break;
	case  1 : return     2*(exponent+1)*(1.0/2.0+(-1.0/2.0+r2/6)*r2); break;
	case  2 : return     2*(exponent+1)*(1.0/2.0+(-3.0/4.0+(1.0/2.0-r2/8)*r2)*r2); break;
	case  3 : return     2*(exponent+1)*(1.0/2.0+(-1.0+(1.0+(-1.0/2.0+r2/10)*r2)*r2)*r2); break;
	case  4 : return     2*(exponent+1)*(1.0/2.0+(-5.0/4.0+(5.0/3.0+(-5.0/4.0+(1.0/2.0-r2/12)*r2)*r2)*r2)*r2); break;
	case  5 : return     2*(exponent+1)*(1.0/2.0+(-3.0/2.0+(5.0/2.0+(-5.0/2.0+(3.0/2.0+(-1.0/2.0+r2/14)*r2)*r2)*r2)*r2)*r2); break;
	case  6 : return     2*(exponent+1)*(1.0/2.0+(-7.0/4.0+(7.0/2.0+(-35.0/8.0+(7.0/2.0+(-7.0/4.0+(1.0/2.0-r2/16)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  7 : return     2*(exponent+1)*(1.0/2.0+(-2.0+(14.0/3.0+(-7.0+(7.0+(-14.0/3.0+(2.0+(-1.0/2.0+r2/18)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  8 : return     2*(exponent+1)*(1.0/2.0+(-9.0/4.0+(6.0+(-21.0/2.0+(63.0/5.0+(-21.0/2.0+(6.0+(-9.0/4.0+(1.0/2.0-r2/20)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  9 : return     2*(exponent+1)*(1.0/2.0+(-5.0/2.0+(15.0/2.0+(-15.0+(21.0+(-21.0+(15.0+(-15.0/2.0+(5.0/2.0+(-1.0/2.0+r2/22)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	}
    return 0.0;
	}

   // Returns 2*Pi*2*q'(r^2) where u(r) = q(r^2)

    double evaluate2xDq2D_4thOrder(double r2) const
	{
	switch(this->exponent-1)
	{
    case  0 : return(-12.0)*(-0.1e1/0.2e1+(0.3e1/0.4e1-r2/0.3e1)*r2);break;
    case  1 : return(-24.0)*(-0.1e1/0.3e1+(0.3e1/0.4e1+(-0.2e1/0.3e1+0.5e1/0.24e2*r2)*r2)*r2);break;
    case  2 : return(-40.0)*(-0.1e1/0.4e1+(0.3e1/0.4e1+(-0.1e1+(0.5e1/0.8e1-0.3e1/0.20e2*r2)*r2)*r2)*r2);break;
	case  3 : return(-60.0)*(-0.1e1/0.5e1+(0.3e1/0.4e1+(-0.4e1/0.3e1+(0.5e1/0.4e1+(-0.3e1/0.5e1+0.7e1/0.60e2*r2)*r2)*r2)*r2)*r2);break;
	case  4 : return(-84.0)*(-0.1e1/0.6e1+(0.3e1/0.4e1+(-0.5e1/0.3e1+(0.25e2/0.12e2+(-0.3e1/0.2e1+(0.7e1/0.12e2-0.2e1/0.21e2*r2)*r2)*r2)*r2)*r2)*r2);break;
	case  5 : return(-112.0)*(-0.1e1/0.7e1+(0.3e1/0.4e1+(-0.2e1+(0.25e2/0.8e1+(-0.3e1+(0.7e1/0.4e1+(-0.4e1/0.7e1+0.9e1/0.112e3*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	case  6 : return(-144.0)*(-0.1e1/0.8e1+(0.3e1/0.4e1+(-0.7e1/0.3e1+(0.35e2/0.8e1+(-0.21e2/0.4e1+(0.49e2/0.12e2+(-0.2e1+(0.9e1/0.16e2-0.5e1/0.72e2*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	case  7 : return(-180.0)*(-0.1e1/0.9e1+(0.3e1/0.4e1+(-0.8e1/0.3e1+(0.35e2/0.6e1+(-0.42e2/0.5e1+(0.49e2/0.6e1+(-0.16e2/0.3e1+(0.9e1/0.4e1+(-0.5e1/0.9e1+0.11e2/0.180e3*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	case  8 : return(-220.0)*(-0.1e1/0.10e2+(0.3e1/0.4e1+(-0.3e1+(0.15e2/0.2e1+(-0.63e2/0.5e1+(0.147e3/0.10e2+(-0.12e2+(0.27e2/0.4e1+(-0.5e1/0.2e1+(0.11e2/0.20e2-0.3e1/0.55e2*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	case  9 : return(-264.0)*(-0.1e1/0.11e2+(0.3e1/0.4e1+(-0.10e2/0.3e1+(0.75e2/0.8e1+(-0.18e2+(0.49e2/0.2e1+(-0.24e2+(0.135e3/0.8e1+(-0.25e2/0.3e1+(0.11e2/0.4e1+(-0.6e1/0.11e2+0.13e2/0.264e3*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	}
    return 0.0;
	}

   // Returns 2*Pi*2*q'(r^2) where u(r) = q(r^2)

    double evaluate2xDq2D_6thOrder(double r2) const
	{
	switch(this->exponent-1)
	{
	case 0 : return(24.0)*(0.1e1/0.2e1+(-0.3e1/0.2e1+(0.5e1/0.3e1-0.5e1/0.8e1*r2)*r2)*r2);break;
    case 1 : return(60.0)*(0.1e1/0.4e1+(-0.1e1+(0.5e1/0.3e1+(-0.5e1/0.4e1+0.7e1/0.20e2*r2)*r2)*r2)*r2);break;
    case 2 : return(120.0)*(0.3e1/0.20e2+(-0.3e1/0.4e1+(0.5e1/0.3e1+(-0.15e2/0.8e1+(0.21e2/0.20e2-0.7e1/0.30e2*r2)*r2)*r2)*r2)*r2);break;
    case 3 : return(210.0)*(0.1e1/0.10e2+(-0.3e1/0.5e1+(0.5e1/0.3e1+(-0.5e1/0.2e1+(0.21e2/0.10e2+(-0.14e2/0.15e2+0.6e1/0.35e2*r2)*r2)*r2)*r2)*r2)*r2);break;
    case 4 : return(336.0)*(0.1e1/0.14e2+(-0.1e1/0.2e1+(0.5e1/0.3e1+(-0.25e2/0.8e1+(0.7e1/0.2e1+(-0.7e1/0.3e1+(0.6e1/0.7e1-0.15e2/0.112e3*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    case 5 : return(504.0)*(0.3e1/0.56e2+(-0.3e1/0.7e1+(0.5e1/0.3e1+(-0.15e2/0.4e1+(0.21e2/0.4e1+(-0.14e2/0.3e1+(0.18e2/0.7e1+(-0.45e2/0.56e2+0.55e2/0.504e3*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    case 6 : return(720.0)*(0.1e1/0.24e2+(-0.3e1/0.8e1+(0.5e1/0.3e1+(-0.35e2/0.8e1+(0.147e3/0.20e2+(-0.49e2/0.6e1+(0.6e1+(-0.45e2/0.16e2+(0.55e2/0.72e2-0.11e2/0.120e3*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    case 7 : return(990.0)*(0.1e1/0.30e2+(-0.1e1/0.3e1+(0.5e1/0.3e1+(-0.5e1+(0.49e2/0.5e1+(-0.196e3/0.15e2+(0.12e2+(-0.15e2/0.2e1+(0.55e2/0.18e2+(-0.11e2/0.15e2+0.13e2/0.165e3*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    case 8 : return(1320.0)*(0.3e1/0.110e3+(-0.3e1/0.10e2+(0.5e1/0.3e1+(-0.45e2/0.8e1+(0.63e2/0.5e1+(-0.98e2/0.5e1+(0.108e3/0.5e1+(-0.135e3/0.8e1+(0.55e2/0.6e1+(-0.33e2/0.10e2+(0.39e2/0.55e2-0.91e2/0.1320e4*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    case 9 : return(1716.0)*(0.1e1/0.44e2+(-0.3e1/0.11e2+(0.5e1/0.3e1+(-0.25e2/0.4e1+(0.63e2/0.4e1+(-0.28e2+(0.36e2+(-0.135e3/0.4e1+(0.275e3/0.12e2+(-0.11e2+(0.39e2/0.11e2+(-0.91e2/0.132e3+0.35e2/0.572e3*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	}
    return 0.0;
	}

    // Creates 2*Pi* value of the derivative for r2 <= 1

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative

    // 0 <=  differentiability order <= 9

    void evaluateDerivatives2D(double r2, double x, double y, vector <double>& derivativeList) const
    {
    	int size = derivativeList.size();
        if(size < 1) return;

        switch(order)
    	{
			case 2  : derivativeList[0] = evaluation2D_2ndOrder(r2); break;
			case 4  : derivativeList[0] = evaluation2D_4thOrder(r2); break;
			case 6  : derivativeList[0] = evaluation2D_6thOrder(r2); break;
    	}
        if(size == 1) {return;}

        double  dVal;

    	switch(order)
    	{
    	case 2 : dVal = evaluate2xDq2D_2ndOrder(r2); break;
    	case 4 : dVal = evaluate2xDq2D_4thOrder(r2); break;
    	case 6 : dVal = evaluate2xDq2D_6thOrder(r2); break;
    	}
    	derivativeList[1] =  dVal*x;
    	derivativeList[2] =  dVal*y;
    }





    double        xPos;    // The position of the mollifier
    double        yPos;

    double         radius;    // The radius of the mollifier
    double   laplaceCoeff;    // Coefficient of the Laplace operator
    double       strength;    // The strength of the mollifier

    int           order;    // Order of the mollifier = 2, 4 or 6
    int        exponent;    // The exponent of the mollifer
    double    logRadius;
};

#undef _DEFAULT_ORDER_
#undef _DEFAULT_SOURCE_DIFFERENTIABLITY_


#endif /* SMOOTHPolyPotential_H_ */


/* Maple code which was used to generate the high order ( > 2) coefficients */
/*
> readlib(C);with(linalg):
> rMin := 0; rMax :=1; iStart:= 4; iMax := 10;
> dim := 3; for M from iStart by 1 to iMax do printf(`case %d : \n`,M-1); qSize:= 2: for q from 0 by 1 to qSize-1 do p := r->a*((1-r^2)^M)+b*((1-r^2)^(M+1));int(r^(dim-1)*r^(2*q)*p(r),r=0..1);v[q+1] := coeffs(int(r^(dim-1)*r^(2*q)*p(r),r=rMin..rMax)); od:w := (matrix(qSize,qSize,[[v[1]],[v[2]]])):iW := inverse(w):cW := multiply(diag(1/iW[1,1],1/iW[1,1]),iW):scaleFactor := iW[1,1]:mcoeff:=array(1..2,[]): mcoeff[1] := cW[1,1]: mcoeff[2] := cW[2,1]:  zzz:=r->(mcoeff[1]*((1-r^2)^M) + mcoeff[2]*((1-r^2)^(M+1))):pr(zzz(r));pr(expand(zzz(r)));  co:=array(1..M+2,[]); co[1] := subs(r=0,zzz(r)); for k from 2 by 1 to M+2 do co[k] := subs(r=0,diff(zzz(r),r$(2*(k-1)))/(2*(k-1))!); od; pr(co); uo:=array(1..M+2,[]);for j from 1 to M+2 do uo[j] := co[j]/((2*j)*(2*j + (dim-2))); od:  pr(uo); fA[1] := uo[1]*r^2: for k from 2 by 1 to M+2 do fA[k] := fA[k-1]+ uo[k]*r^(2*k): od: u := r->convert(fA[M+2],polynom):pr(u(r)-subs(r=1,u(r))); pr(simplify((1/((r^(dim-1)))*diff(r^(dim-1)*diff(u(r),r),r)))); pr(expand(zzz(r))); print(C(scaleFactor));print(C(convert(u(r)-subs(r=1,u(r)),horner),optimized)); od:

> qSize:= 3: for M from iStart by 1 to iMax do printf(`case %d : \n`,M-1); for q from 0 by 1 to qSize-1 do p := r->a*((1-r^2)^M)+b*((1-r^2)^(M+1)) + c*((1-r^2)^(M+2));int(r^(dim-1)*r^(2*q)*p(r),r=rMin..rMax);v[q+1] := coeffs(int(r^(dim-1)*r^(2*q)*p(r),r=rMin..rMax)); od:w := (matrix(qSize,qSize,[[v[1]],[v[2]],[v[3]]])):iW := inverse(w):cW := multiply(diag(1/iW[1,1],1/iW[1,1],1/iW[1,1]),iW):scaleFactor := iW[1,1]:mcoeff:=array(1..qSize,[]): mcoeff[1] := cW[1,1]: mcoeff[2] := cW[2,1]:mcoeff[3] := cW[3,1]: zzz:=r->(mcoeff[1]*((1-r^2)^M) + mcoeff[2]*((1-r^2)^(M+1)) + mcoeff[3]*((1-r^2)^(M+2))):pr(zzz(r));pr(expand(zzz(r)));  co:=array(1..M+3,[]); co[1] := subs(r=0,zzz(r)); for k from 2 by 1 to M+3 do co[k] := subs(r=0,diff(zzz(r),r$(2*(k-1)))/(2*(k-1))!); od; pr(co); uo:=array(1..M+3,[]);for j from 1 to M+3 do uo[j] := co[j]/((2*j)*(2*j + (dim-2))); od:  pr(uo); fA[1] := (uo[1])*r^2: for k from 2 by 1 to M+3 do fA[k] := fA[k-1]+ uo[k]*r^(2*k): od: u := r->convert(fA[M+3],polynom):pr(u(r)); pr(simplify((1/((r^(dim-1)))*diff(r^(dim-1)*diff(u(r),r),r)))); pr(expand(zzz(r))); print(C(scaleFactor)); print(C(convert(u(r) - subs(r=1,u(r)),horner),optimized)); od:


*/
