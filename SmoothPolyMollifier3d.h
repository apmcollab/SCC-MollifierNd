/*
 * SmoothPolyMollifier3d.h
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
#include <vector>
#include <functional>

#ifndef SMOOTH_POLY_MOLLIFIER_3D_
#define SMOOTH_POLY_MOLLIFIER_3D_

#define DEFAULT_DIFFERENTIABLITY 6

class SmoothPolyMollifier3d
{
	public:

	SmoothPolyMollifier3d()
	{
	initialize();
	}

    SmoothPolyMollifier3d(const SmoothPolyMollifier3d& S)
	{
    initialize(S);
	}

    SmoothPolyMollifier3d(double xPos, double yPos,double zPos, double radius, double strength)
	{
    initialize(xPos,yPos,zPos,radius,strength);
	}

	void initialize()
	{
    radius     = 0.0;
    strength   = 0.0;
    xPos       = 0.0;
    yPos       = 0.0;
    zPos       = 0.0;
    exponent   = DEFAULT_DIFFERENTIABLITY + 1 ;
	}

	void initialize(const SmoothPolyMollifier3d& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    yPos       = S.yPos;
    zPos       = S.zPos;
    radius     = S.radius;
    exponent   = S.exponent;
	}


	void initialize(double xPos, double yPos, double zPos, double radius, double strength)
	{
    this->exponent   = DEFAULT_DIFFERENTIABLITY + 1;

    this->strength   = strength;
    this->radius     = radius;

    this->xPos       = xPos;
    this->yPos       = yPos;
    this->zPos       = zPos;
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


	double operator()(double x, double y, double z) const
	{
    double r2radius = ((x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos))/(radius*radius);
    if(r2radius >= 1.0) return 0.0;
    return (strength/(radius*radius*radius))*evaluation3D(r2radius);
	}

    //
    // Normalized evaluation routine:
    //
    // Input:  A value r2 = r*r for r in [0,1]
    //
    // Output:  Returns (sigma/(4*Pi))*(1-r2)^(exponent)
    // where sigma is the prefactor so that the mollifier
    // has unit integral over the unit ball.
    //

    double evaluation3D(double r2) const
	{
    double u = 1.0-r2;
    double oneOver4pi = 7.9577471545947667884e-02;
	switch(this->exponent-1)
    {
	case  0 : return   oneOver4pi*(15.0/2.0)*u;  break;
    case  1 : return   oneOver4pi*(105.0/8.0)*u*u; break;
    case  2 : return   oneOver4pi*(315.0/16.0)*u*u*u; break;
    case  3 : return   oneOver4pi*(3465.0/128.0)*u*u*u*u; break;
    case  4 : return   oneOver4pi*(9009.0/256.0)*u*u*u*u*u; break;
    case  5 : return   oneOver4pi*(45045.0/1024.0)*u*u*u*u*u*u; break;
    case  6 : return   oneOver4pi*(109395.0/2048.0)*u*u*u*u*u*u*u; break;
    case  7 : return   oneOver4pi*(2078505.0/32768.0)*u*u*u*u*u*u*u*u; break;
    case  8 : return   oneOver4pi*(4849845.0/65536.0)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return   oneOver4pi*(22309287.0/262144.0)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
	}


	 //  Returns a std::function that is bound to the evaluation operator of *this

	std::function<double(double,double,double)> getEvaluationPtr() const
	{
	std::function<double(double,double,double)> F = [this](double x,double y,double z) {return this->operator()(x,y,z);};
	return F;
	}

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
    double r2radius = xRadius*xRadius + yRadius*yRadius + zRadius*zRadius;

    double oneOver4pi = 7.9577471545947667884e-02;

    int dCount = (int)derivativeList.size();

    if(r2radius >= 1.0)
    {
    for(long i = 0; i < dCount; i++)
    {
    derivativeList[i] = 0.0;
    }
    return;
    }

    evaluateDerivatives3D(r2radius, xRadius, yRadius, zRadius, derivativeList);

    if(dCount >= 1)
    {
    derivativeList[0] *= oneOver4pi*(strength/(radius*radius*radius));
    }

    if(dCount >= 4)
    {
    	for(long i = 1; i <= 3; i++)
    	{
    	derivativeList[i] *= oneOver4pi*(strength/(radius*radius*radius*radius));
    	}
    }

    if(dCount >= 10)
    {
    	for(long i = 4; i <= 9; i++)
    	{
    	derivativeList[i] *= oneOver4pi*(strength/(radius*radius*radius*radius*radius));
    	}
    }
    return;
   	}


private:

	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	//    Utility functions for evaluating the derivatives of the mollifer
	//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    // Creates 4*Pi* value of the derivative of for r2 <= 1

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

    // u(r)         = sigma*(1-r^2)^M
	// du/dx(r)     = sigma*M*(1-r^2)^(M-1)*(-2 r) (x/r) = -2.0*sigma*M*(1-r^2)^(M-1)*x
	// d^2u/dx^2(r) = 4.0*sigma*M*(M-1)*(1-r^2)^(M-2)*x*x - 2.0*sigma*M*(1-r^2)^(M-1)

    void evaluateDerivatives3D(double r2, double x, double y, double z, std::vector <double>& derivativeList) const
    {
    	int size = (int)derivativeList.size();
    	if(size == 1)  {derivativeList[0] = 1.2566370614359172954e+01*evaluation3D(r2); return;}

    	double sigma = 0.0;

		switch(this->exponent-1)
		{
		case  0 : sigma =     (15.0/2.0); break;
		case  1 : sigma =     (105.0/8.0);  break;
		case  2 : sigma =     (315.0/16.0); break;
		case  3 : sigma =     (3465.0/128.0); break;
		case  4 : sigma =     (9009.0/256.0); break;
		case  5 : sigma =     (45045.0/1024.0); break;
		case  6 : sigma =     (109395.0/2048.0); break;
		case  7 : sigma =     (2078505.0/32768.0); break;
		case  8 : sigma =     (4849845.0/65536.0); break;
		case  9 : sigma =     (22309287.0/262144.0); break;
		}

		double u     = 1.0-r2;

    	double uMM1  = 1.0; // (1-r^2)^(exponent-1)
    	double uMM2  = 1.0; // (1-r^2)^(exponent-2)
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
    	if(exponent-2 == -1)
    	{
    	uMM2 = 0.0;
    	uMM1 = 1.0;
    	}

    	if(size >= 4)
    	{
    	derivativeList[0] =  sigma*uMM1*u;
    	derivativeList[1] = -2.0*sigma*uMM1*M*x;
    	derivativeList[2] = -2.0*sigma*uMM1*M*y;
    	derivativeList[3] = -2.0*sigma*uMM1*M*z;
    	}

    	double MM2factor = 4.0*sigma*M*(M-1.0);
    	if(size >= 10)
    	{
    	derivativeList[4] = MM2factor*uMM2*x*x - 2.0*sigma*uMM1*M;
    	derivativeList[5] = MM2factor*uMM2*x*y;
    	derivativeList[6] = MM2factor*uMM2*x*z;
    	derivativeList[7] = MM2factor*uMM2*y*y - 2.0*sigma*uMM1*M;
    	derivativeList[8] = MM2factor*uMM2*y*z;
    	derivativeList[9] = MM2factor*uMM2*z*z - 2.0*sigma*uMM1*M;
    	}
    }


    double    strength;    // The strength of the mollifier
    double        xPos;    // The position of the mollifier
    double        yPos;
    double        zPos;

    double        radius;   // The radius of the mollifier
    int         exponent;   // The exponent of the mollifier
};


#undef DEFAULT_DIFFERENTIABLITY
#endif /* SMOOTHPOLYMOLLIFIER_H_ */
