/*
 * SmoothPolyPotential3d.h
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

#include "SmoothPolyMollifier3d.h"

#ifndef SMOOTH_POLY_POTENTIAL_3D_
#define SMOOTH_POLY_POTENTIAL_3D_

#define _DEFAULT_SOURCE_DIFFERENTIABLITY_ 6

/*
 	A class whose member functions

 	operator(double x, double y, double z) (3D)

 	evaluate the potential based upon source mollifiers of the form

 	(SIGMA_M(N)/OMEGA_N) ( 1 - r^2)^exponent

 	for exponent = 1..10 and centered at

 	(xPos,yPos,zPos) of specified strength and radius.

 	This class also contains routines for evaluating the
 	first and second derivatives of the potential.

 */
class SmoothPolyPotential3d
{
	public:

	SmoothPolyPotential3d()
	{
	initialize();
	}

	SmoothPolyPotential3d(const SmoothPolyPotential3d& S)
	{
    initialize(S);
	}


    SmoothPolyPotential3d(double xPos, double yPos,double zPos, double radius, double strength, double laplaceCoeff = 1.0)
	{
    initialize(xPos, yPos, zPos ,radius, strength,laplaceCoeff);
	}

	void initialize()
	{
	xPos       = 0.0;
	yPos       = 0.0;
	zPos       = 0.0;
    radius     = 0.0;
    laplaceCoeff  = 1.0;
    strength      = 0.0;
    exponent      = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;

    source.initialize();
	}

    void initialize(const SmoothPolyPotential3d& S)
	{
    xPos       = S.xPos;
    yPos       = S.yPos;
    zPos       = S.zPos;
    radius     = S.radius;
    laplaceCoeff   = S.laplaceCoeff;
    strength       = S.strength;
    exponent       = S.exponent;

    source.initialize(S.source);
	}

	void initialize(double xPos, double yPos, double zPos, double radius, double strength, double laplaceCoeff = 1.0)
	{
	this->xPos      = xPos;
    this->yPos      = yPos;
    this->zPos      = zPos;
    this->radius    = radius;
    this->laplaceCoeff = laplaceCoeff;
    this->strength     = strength;
    this->exponent     = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;

    source.initialize(xPos,yPos,zPos,radius,strength);
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

	SmoothPolyMollifier3d getSource() const
	{
	return source;
	}

    void setLaplaceCoefficient(double laplaceCoeff)
	{
	this->laplaceCoeff = laplaceCoeff;
	}

	double getLaplaceCoefficient() const
	{
	return this->laplaceCoeff;
	}

    double operator()(double x, double y, double z) const
	{
	double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos);
	double r2radius = r2/(radius*radius);
	if(r2radius >= 1.0)  {return -(strength/laplaceCoeff)*0.79577471545948e-01*(1.0/(std::sqrt(r2)));} // strength*(-1.0/(4*pi*r)
    return ((evaluation3D_2ndOrder(r2radius) - 1.0)*.79577471545947667883e-01*(strength/laplaceCoeff))/radius;
	}

	//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double,double,double )> getEvaluationPtr() const
	{
	std::function<double(double,double,double)> F = [this](double x,double y, double z) {return this->operator()(x,y,z);};
	return F;
	}



	double evaluateSource(double x, double y, double z) const
	{
	return source.operator()(x,y,z);
	}


	std::function<double(double,double,double )> getSourceEvaluationPtr() const
	{
	std::function<double(double,double,double)> F = [this](double x,double y, double z) {return this->evaluateSource(x,y,z);};
	return F;
	}


    // returns 4*pi* value

    double evaluation3D_2ndOrder(double r2) const
	{
	switch(this->exponent-1)
	{
	case  0 : return    (15.0/2.0)*((1.0/6.0-r2/20)*r2 -  7.0/60.0); break;
	case  1 : return    (105.0/8.0)*((1.0/6.0+(-1.0/10.0+r2/42)*r2)*r2 -19.0/210.0);  break;
	case  2 : return    (315.0/16.0)*((1.0/6.0+(-3.0/20.0+(1.0/14.0-r2/72.0)*r2)*r2)*r2 - 187.0/2520.0); break;
	case  3 : return    (3465.0/128.0)*((1.0/6.0+(-1.0/5.0+(1.0/7.0+(-1.0/18.0+r2/110.0)*r2)*r2)*r2)*r2 - 437.0/6930.0); break;
	case  4 : return    (9009.0/256.0)*((1.0/6.0+(-1.0/4.0+(5.0/21.0+(-5.0/36.0+(1.0/22.0-r2/156.0)*r2)*r2)*r2)*r2)*r2 - 1979.0/36036.0); break;
	case  5 : return    (45045.0/1024.0)*((1.0/6.0+(-3.0/10.0+(5.0/14.0+(-5.0/18.0+(3.0/22.0+(-1.0/26.0+r2/210)*r2)*r2)*r2)*r2)*r2)*r2 - 4387.0/90090.0); break;
	case  6 : return    (109395.0/2048.0)*((1.0/6.0+(-7.0/20.0+(1.0/2.0+(-35.0/72.0+(7.0/22.0+(-7.0/52.0+(1.0/30.0-r2/272.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 76627.0/1750320.0); break;
	case  7 : return    (2078505.0/32768.0)*((1.0/6.0+(-2.0/5.0+(2.0/3.0+(-7.0/9.0+(7.0/11.0+(-14.0/39.0+(2.0/15.0+(-1.0/34.0+r2/342.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 165409.0/4157010.0); break;
	case  8 : return    (4849845.0/65536.0)*((1.0/6.0+(-9.0/20.0+(6.0/7.0+(-7.0/6.0+(63.0/55.0+(-21.0/26.0+(2.0/5.0+(-9.0/68.0+(1.0/38.0-r2/420)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 141565.0/3879876.0); break;
	case  9 : return    (22309287.0/262144.0)*((1.0/6.0+(-1.0/2.0+(15.0/14.0+(-5.0/3.0+(21.0/11.0+(-21.0/13.0+(1.0+(-15.0/34.0+(5.0/38.0+(-1.0/42.0+r2/506.0)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2 - 1503829.0/44618574.0); break;
	}
    return 0.0;
	}

    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =   z derivative
    // derivativeList[4] =  xx derivative
    // derivativeList[5] =  xy derivative
    // derivativeList[6] =  xz derivative
    // derivativeList[7] =  yy derivative
    // derivativeList[8] =  yz derivative
    // derivativeList[9] =  zz derivative
    //

   	void derivatives(double x, double y, double z, std::vector <double>& derivativeList, int maxOrder = 2) const
   	{

   	switch (maxOrder)
   	{
   	case 0  :  derivativeList.resize(1); break;
   	case 1  :  derivativeList.resize(4); break;
   	case 2  :  derivativeList.resize(10); break;
   	default :  derivativeList.resize(10); break;
   	}

   	double xRadius  = (x-xPos)/radius;
   	double yRadius  = (y-yPos)/radius;
    double zRadius  = (z-zPos)/radius;
    double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos) ;
	double r2radius = r2/(radius*radius);
	double r;
	double r3;
	double r5;

    int dCount = (int)derivativeList.size();

    if(r2radius <= 1.0)
    {
    evaluateDerivatives3D(r2radius, xRadius, yRadius, zRadius,derivativeList);
    }
    else
    {
    	r = std::sqrt(r2);
    	if(dCount >= 1)
    	{
    		derivativeList[0] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*(1.0/r);
    	}
    	if(dCount >= 4)
    	{
    	    r3 = r2*r;
    		derivativeList[1] =  (strength/laplaceCoeff)*.79577471545947667883e-01*((x-xPos)/(r3));
    		derivativeList[2] =  (strength/laplaceCoeff)*.79577471545947667883e-01*((y-yPos)/(r3));
    		derivativeList[3] =  (strength/laplaceCoeff)*.79577471545947667883e-01*((z-zPos)/(r3));
    	}
    	if(dCount >= 10)
    	{
    	    r5 = r3*r2;
    		derivativeList[4] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*((3.0*(x-xPos)*(x-xPos))/(r5) - (1.0/(r3)));
    		derivativeList[5] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*((3.0*(x-xPos)*(y-yPos))/(r5));
    		derivativeList[6] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*((3.0*(x-xPos)*(z-zPos))/(r5));
    		derivativeList[7] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*((3.0*(y-yPos)*(y-yPos))/(r5) - (1.0/(r3)));
    		derivativeList[8] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*((3.0*(y-yPos)*(z-zPos))/(r5));
    		derivativeList[9] =  -(strength/laplaceCoeff)*.79577471545947667883e-01*((3.0*(z-zPos)*(z-zPos))/(r5) - (1.0/(r3)));
    	}
    return;
    }

    // Scale results

    if(dCount >= 1)
    {
    derivativeList[0] -= 1.0;
    derivativeList[0] *= ((strength/laplaceCoeff)*.79577471545947667883e-01)/radius;
    }

    if(dCount >= 4)
    {
    	for(long i = 1; i <= 3; i++)
    	{
    	derivativeList[i] *= ((strength/laplaceCoeff)*.79577471545947667883e-01)/(radius*radius);
    	}
    }

    if(dCount >= 10)
    {
    	for(long i = 4; i <= 9; i++)
    	{
    	derivativeList[i] *= ((strength/laplaceCoeff)*.79577471545947667883e-01)/(radius*radius*radius);
    	}
    }

    return;
   	}


   	double radialDerivative(double x, double y, double z) const
   	{
    double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos) ;
	double r2radius = r2/(radius*radius);
	double r = std::sqrt(r2);

	if(r2radius > 1.0)
	{
		return  ((strength/laplaceCoeff)*.79577471545947667883e-01)/r2;
	}

	double dVal = evaluate2xDq3D(r2radius);

	dVal *=  (r*(strength/laplaceCoeff)*.79577471545947667883e-01)/(radius*radius*radius);
	return dVal;
   	}

	private:
	//
	// Routines required for the evaluation of the derivatives of the potentials
	//

    // Returns 4*Pi*2*q'(r^2) where u(r) = q(r^2)

    double evaluate2xDq3D(double r2) const
	{
	switch(this->exponent-1)
	{
	case  0 : return    (15.0/2.0)*(1.0/3.0-r2/5); break;
	case  1 : return    (105.0/8.0)*(1.0/3.0+(-2.0/5.0+r2/7)*r2);  break;
	case  2 : return    (315.0/16.0)*(1.0/3.0+(-3.0/5.0+(3.0/7.0-r2/9)*r2)*r2); break;
	case  3 : return    (3465.0/128.0)*(1.0/3.0+(-4.0/5.0+(6.0/7.0+(-4.0/9.0+r2/11)*r2)*r2)*r2); break;
	case  4 : return    (9009.0/256.0)*(1.0/3.0+(-1.0+(10.0/7.0+(-10.0/9.0+(5.0/11.0-r2/13)*r2)*r2)*r2)*r2); break;
	case  5 : return    (45045.0/1024.0)*(1.0/3.0+(-6.0/5.0+(15.0/7.0+(-20.0/9.0+(15.0/11.0+(-6.0/13.0+r2/15)*r2)*r2)*r2)*r2)*r2); break;
	case  6 : return    (109395.0/2048.0)*(1.0/3.0+(-7.0/5.0+(3.0+(-35.0/9.0+(35.0/11.0+(-21.0/13.0+(7.0/15.0-r2/17)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  7 : return    (2078505.0/32768.0)*(1.0/3.0+(-8.0/5.0+(4.0+(-56.0/9.0+(70.0/11.0+(-56.0/13.0+(28.0/15.0+(-8.0/17.0+r2/19)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  8 : return    (4849845.0/65536.0)*(1.0/3.0+(-9.0/5.0+(36.0/7.0+(-28.0/3.0+(126.0/11.0+(-126.0/13.0+(28.0/5.0+(-36.0/17.0+(9.0/19.0-r2/21)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  9 : return    (22309287.0/262144.0)*(1.0/3.0+(-2.0+(45.0/7.0+(-40.0/3.0+(210.0/11.0+(-252.0/13.0+(14.0+(-120.0/17.0+(45.0/19.0+(-10.0/21.0+r2/23)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	}
    return 0.0;
	}

	//  Returns 4*Pi*(4*q''(r^2)) where u(r) = q(r^2)

	double evaluate4xD2q3D(double r2) const
	{
	switch(this->exponent-1)
	{
	case  0 : return    (15.0/2.0)*(-2.0/5.0); break;
	case  1 : return    (105.0/8.0)*(-4.0/5.0+4.0/7.0*r2);  break;
	case  2 : return    (315.0/16.0)*(-6.0/5.0+(12.0/7.0-2.0/3.0*r2)*r2); break;
	case  3 : return    (3465.0/128.0)*(-8.0/5.0+(24.0/7.0+(-8.0/3.0+8.0/11.0*r2)*r2)*r2); break;
	case  4 : return    (9009.0/256.0)*(-2.0+(40.0/7.0+(-20.0/3.0+(40.0/11.0-10.0/13.0*r2)*r2)*r2)*r2); break;
	case  5 : return    (45045.0/1024.0)*(-12.0/5.0+(60.0/7.0+(-40.0/3.0+(120.0/11.0+(-60.0/13.0+4.0/5.0*r2)*r2)*r2)*r2)*r2); break;
	case  6 : return    (109395.0/2048.0)*( -14.0/5.0+(12.0+(-70.0/3.0+(280.0/11.0+(-210.0/13.0+(28.0/5.0-14.0/17.0*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  7 : return    (2078505.0/32768.0)*(-16.0/5.0+(16.0+(-112.0/3.0+(560.0/11.0+(-560.0/13.0+(112.0/5.0+(-112.0/17.0+16.0/19.0*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  8 : return    (4849845.0/65536.0)*(-18.0/5.0+(144.0/7.0+(-56.0+(1008.0/11.0+(-1260.0/13.0+(336.0/5.0+(-504.0/17.0+(144.0/19.0-6.0/7.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	case  9 : return    (22309287.0/262144.0)*(-4.0+(180.0/7.0+(-80.0+(1680.0/11.0+(-2520.0/13.0+(168.0+(-1680.0/17.0+(720.0/19.0+(-60.0/7.0+20.0/23.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
	}
    return 0.0;
	}


    // Creates 4*Pi* value of the derivative for r2 <= 1

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =   z derivative
    // derivativeList[4] =  xx derivative
    // derivativeList[5] =  xy derivative
    // derivativeList[6] =  xz derivative
    // derivativeList[7] =  yy derivative
    // derivativeList[8] =  yz derivative
    // derivativeList[9] =  zz derivative
    //
    // 0 <=  differentiability order <= 9

    void evaluateDerivatives3D(double r2, double x, double y, double z, std::vector <double>& derivativeList) const
    {
    	int size = (int)derivativeList.size();
    	if(size == 1)  {derivativeList[0] = evaluation3D_2ndOrder(r2); return;}

        double  dVal;
        double d2Val;
    	if(size >= 4)
    	{
    	dVal = evaluate2xDq3D(r2);
    	derivativeList[0] =  evaluation3D_2ndOrder(r2);
    	derivativeList[1] =  dVal*x;
    	derivativeList[2] =  dVal*y;
    	derivativeList[3] =  dVal*z;
    	}

    	if(size >= 10)
    	{
    	d2Val = evaluate4xD2q3D(r2);
    	derivativeList[4] = d2Val*x*x + dVal;
    	derivativeList[5] = d2Val*x*y;
    	derivativeList[6] = d2Val*x*z;
    	derivativeList[7] = d2Val*y*y + dVal;
    	derivativeList[8] = d2Val*y*z;
    	derivativeList[9] = d2Val*z*z + dVal;
    	}
    }


    double    strength;    // The strength of the mollifier
    double        xPos;    // The position of the mollifier
    double        yPos;
    double        zPos;

    double  laplaceCoeff;     // Coefficient of the Laplacian
    double        radius;   // The radius of the mollifier
    int         exponent;   // The exponent of the mollifier, e.g. (1-r^2)^exponent.

    SmoothPolyMollifier3d source;

};

#undef _DEFAULT_SOURCE_DIFFERENTIABLITY_


#endif /* SMOOTHPolyPotential_H_ */
