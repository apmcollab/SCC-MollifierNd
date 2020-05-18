/*
 * SmoothPolyPotential2d.h
 *
 *  Created on: Feb 16, 2014
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

#include "SmoothPolyMollifier2d.h"

#ifndef _SmoothPolyPotential2d_
#define _SmoothPolyPotential2d_

#define _DEFAULT_SOURCE_DIFFERENTIABLITY_ 6

/*
 	A class whose member function

 	operator(double x, double y)

 	evaluate the potential based upon mollifiers of the form

 	(SIGMA_M(2)/OMEGA_2) ( 1 - r^2)^exponent

 	for exponent = 1..10 and centered at

 	(xPos,yPos) of specified strength and radius.

 	This class also contains routines for evaluating the
 	first and second derivatives of the potential.

 */
class SmoothPolyPotential2d
{
	public:

	SmoothPolyPotential2d()
	{
	initialize();
	}

	SmoothPolyPotential2d(const SmoothPolyPotential2d& S)
	{
    initialize(S);
	}

    SmoothPolyPotential2d(double xPos, double yPos, double radius, double strength, double laplaceCoeff = 1.0)
	{
    initialize(xPos,yPos,radius,strength, laplaceCoeff);
	}

	void initialize()
	{
	xPos      = 0.0;
	yPos      = 0.0;
    radius    = 0.0;
    strength  = 0.0;
    exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
    logRadius = 0.0;
    laplaceCoeff      = 1.0;

    source.initialize();
	}

    void initialize(const SmoothPolyPotential2d& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    yPos       = S.yPos;
    radius     = S.radius;
    exponent   = S.exponent;
    logRadius  = S.logRadius;
    laplaceCoeff       = S.laplaceCoeff;
    source.initialize(S.source);
	}


	void initialize(double xPos, double yPos,double radius, double strength, double laplaceCoeff = 1.0)
	{
	this->xPos      = xPos;
    this->yPos      = yPos;
    this->radius    = radius;
    this->strength  = strength;
    this->laplaceCoeff      = laplaceCoeff;

    this->exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
    this->logRadius = std::log(radius);

    source.initialize(xPos,yPos,radius,strength);
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

	SmoothPolyMollifier2d getSource() const
	{
	return source;
	}

	double operator()(double x, double y) const
	{
	double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos);
	double r2radius = r2/(radius*radius);
	if(r2radius >= 1.0) {return (strength/laplaceCoeff)*.15915494309190*std::log(std::sqrt(r2));}   // strength*(1/(2*pi))*std::log(r)
    return (strength/laplaceCoeff)*.15915494309190*(evaluation2D_2ndOrder(r2radius) + logRadius);
	}
	//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double,double)> getEvaluationPtr() const
	{
    std::function<double(double,double)> F = [this](double x,double y) {return this->operator()(x,y);};
	return F;
	}



	double evaluateSource(double x,double y) const
	{
    return source.operator()(x,y);
	}

	//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double,double)> getSourceEvaluationPtr() const
	{
    std::function<double(double,double)> F = [this](double x, double y) {return this->evaluateSource(x,y);};
	return F;
	}



	void setLaplaceCoefficient(double laplaceCoeff)
	{
	this->laplaceCoeff = laplaceCoeff;
	}

	double getLaplaceCoefficient() const
	{
	return this->laplaceCoeff;
	}


	double evaluation2D_2ndOrder(double r2)  const
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

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =  xx derivative
    // derivativeList[4] =  xy derivative
    // derivativeList[5] =  yy derivative
    //

	void derivatives(double x, double y, std::vector <double>& derivativeList, int maxOrder = 2) const
   	{
   	switch (maxOrder)
   	{
   	case 0  :  derivativeList.resize(1); break;
   	case 1  :  derivativeList.resize(3); break;
   	case 2  :  derivativeList.resize(6); break;
   	default :  derivativeList.resize(6); break;
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
    	r = std::sqrt(r2);
    	if(dCount >= 1)
    	{
    		derivativeList[0] =  (strength/laplaceCoeff)*.15915494309190*std::log(r);
    	}
    	if(dCount >= 3)
    	{
    		derivativeList[1] =  (strength/laplaceCoeff)*.15915494309190*((x-xPos)/r2);
    		derivativeList[2] =  (strength/laplaceCoeff)*.15915494309190*((y-yPos)/r2);
    	}
    	if(dCount >= 6)
    	{
    		derivativeList[3] =  (strength/laplaceCoeff)*.15915494309190*((-2.0*(x-xPos)*(x-xPos))/(r2*r2) + (1.0/r2));
    		derivativeList[4] =  (strength/laplaceCoeff)*.15915494309190*((-2.0*(x-xPos)*(y-yPos))/(r2*r2));
    		derivativeList[5] =  (strength/laplaceCoeff)*.15915494309190*((-2.0*(y-yPos)*(y-yPos))/(r2*r2)+ (1.0/r2));
    	}
    return;
    }

    // Scale results

    if(dCount >= 1)
    {
    derivativeList[0] += logRadius;
    derivativeList[0] *= (strength/laplaceCoeff)*.15915494309190;
    }

    if(dCount >= 3)
    {
    	for(long i = 1; i <= 2; i++)
    	{
    	derivativeList[i] *= (strength)/(laplaceCoeff*radius*6.28318530717958648);
    	}
    }

    if(dCount >= 6)
    {
    	for(long i = 3; i <= 5; i++)
    	{
    	derivativeList[i] *= (strength)/(laplaceCoeff*radius*radius*6.28318530717958648);
    	}
    }

    return;
   	}


   	double radialDerivative(double x, double y) const
   	{
    double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos);
	double r2radius = r2/(radius*radius);
	double r = std::sqrt(r2);

	if(r2radius > 1.0)
	{
		return ((strength/laplaceCoeff)*.15915494309190)/r;
	}

	double dVal = evaluate2xDq2D(r2radius);

	dVal *=  (r*strength)/(laplaceCoeff*radius*radius*6.28318530717958648);

	return dVal;
   	}

    private:

	//
	// Routines required for the evaluation of the derivatives of the potentials
	//

    // Returns 2*Pi*2*q'(r^2) where u(r) = q(r^2)

	double evaluate2xDq2D(double r2) const
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

   // Returns 2*Pi*(4*q''(r^2)) where u(r) = q(r^2)

   double evaluate4xD2q2D(double r2) const
   {
   	switch(this->exponent-1)
	{
	case  0 : return     2*(exponent+1)*(-1.0/2.0); break;
	case  1 : return     2*(exponent+1)*(-1.0+2.0/3.0*r2); break;
	case  2 : return     2*(exponent+1)*(-3.0/2.0+(2.0-3.0/4.0*r2)*r2); break;
	case  3 : return     2*(exponent+1)*(-2.0+(4.0+(-3.0+4.0/5.0*r2)*r2)*r2); break;
	case  4 : return     2*(exponent+1)*(-5.0/2.0+(20.0/3.0+(-15.0/2.0+(4.0-5.0/6.0*r2)*r2)*r2)*r2); break;
	case  5 : return     2*(exponent+1)*(-3.0+(10.0+(-15.0+(12.0+(-5.0+6.0/7.0*r2)*r2)*r2)*r2)*r2);   break;
	case  6 : return     2*(exponent+1)*(-7.0/2.0+(14.0+(-105.0/4.0+(28.0+(-35.0/2.0+(6.0-7.0/8.0*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  7 : return     2*(exponent+1)*(-4.0+(56.0/3.0+(-42.0+(56.0+(-140.0/3.0+(24.0+(-7.0+8.0/9.0*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  8 : return     2*(exponent+1)*(-9.0/2.0+(24.0+(-63.0+(504.0/5.0+(-105.0+(72.0+(-63.0/2.0+(8.0-9.0/10.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  9 : return     2*(exponent+1)*(-5.0+(30.0+(-90.0+(168.0+(-210.0+(180.0+(-105.0+(40.0+(-9.0+10.0/11.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
	}
    return 0.0;
   }

    // Creates 2*Pi* value of the derivative for r2 <= 1

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =  xx derivative
    // derivativeList[4] =  xy derivative
    // derivativeList[5] =  yy derivative
    //
    // 0 <=  differentiability order <= 9

    void evaluateDerivatives2D(double r2, double x, double y, std::vector <double>& derivativeList) const
    {
    	int size = derivativeList.size();
    	if(size == 1)  {derivativeList[0] = evaluation2D_2ndOrder(r2); return;}

        double  dVal;
        double d2Val;
    	if(size >= 3)
    	{
    	dVal = evaluate2xDq2D(r2);
    	derivativeList[0] =  evaluation2D_2ndOrder(r2);
    	derivativeList[1] =  dVal*x;
    	derivativeList[2] =  dVal*y;
    	}

    	if(size >= 6)
    	{
    	d2Val = evaluate4xD2q2D(r2);
    	derivativeList[3] = d2Val*x*x + dVal;
    	derivativeList[4] = d2Val*x*y;
    	derivativeList[5] = d2Val*y*y + dVal;
    	}
    }



    double        xPos;    // The position of the mollifier
    double        yPos;
    double    strength;    // The strength of the mollifier

    double        radius;   // The radius of the mollifier
    int         exponent;   // The exponent of the mollifier, e.g. (1-r^2)^exponent.
    double     logRadius;

    double  laplaceCoeff;   // Coefficient of Laplace operator

    SmoothPolyMollifier2d source;
};

#undef _DEFAULT_SOURCE_DIFFERENTIABLITY_


#endif /* SMOOTHPolyPotential_H_ */
