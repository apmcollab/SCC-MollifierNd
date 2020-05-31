/*
 * HighOrderPolyMollifier3d.h
 *
 *  Created on: Feb 16, 2014
 *      Author: anderson
 *
 *  Updated March 11, 2015
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


#ifndef HIGH_ORDER_POLY_MOLLIFIER_3D_
#define HIGH_ORDER_POLY_MOLLIFIER_3D_

#define MOLLIFIER_DEFAULT_ORDER      4
#define DEFAULT_DIFFERENTIABLITY 6

class HighOrderPolyMollifier3d
{
	public:

	HighOrderPolyMollifier3d()
	{
	initialize();
	}

	HighOrderPolyMollifier3d(const HighOrderPolyMollifier3d& S)
	{
    xPos       = S.xPos;
    yPos       = S.yPos;
    zPos       = S.zPos;
    radius     = S.radius;
    strength   = S.strength;

    order      = S.order;
    exponent   = S.exponent;
	}


    HighOrderPolyMollifier3d(double xPos, double yPos,double zPos, double radius, double strength)
	{
    initialize(xPos,yPos,zPos,radius,strength);
	}

	void initialize()
	{
    xPos       = 0.0;
    yPos       = 0.0;
    zPos       = 0.0;
    radius     = 0.0;
    strength   =  0.0;

    order     = MOLLIFIER_DEFAULT_ORDER;
    exponent  = DEFAULT_DIFFERENTIABLITY + 1;
	}

	void initialize(double xPos, double yPos, double zPos, double radius, double strength)
	{
	this->xPos      = xPos;
    this->yPos      = yPos;
    this->zPos      = zPos;

    this->radius    = radius;
    this->strength  = strength;

    this->exponent  = DEFAULT_DIFFERENTIABLITY + 1;
    this->order     = MOLLIFIER_DEFAULT_ORDER;
	}

	void setRadius(double radius) {this->radius = radius;}
	double getRadius()            {return this->radius;  }

	void setOrder(int order)
	{
	this->order = order;
	if((this->order != 2)&&(this->order != 4)&&(this->order != 6))
	{
	this->order = MOLLIFIER_DEFAULT_ORDER;
	}
	}

	int getOrder() const
	{return this->order;}

	// 0 <= diffOrder <= 9

	void setDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;
	}

	int getDifferentiablity() const
	{
	return this->exponent-1;
	}

	/// Evaluation operator 3D

	double operator()(double x, double y, double z) const
	{
    double r2radius = ((x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos))/(radius*radius);
    if(r2radius >= 1.0) return 0.0;
    switch(order)
    {
    case 2 : return (strength/(radius*radius*radius*12.5663706143591730))*evaluation3D_2ndOrder(r2radius);
    case 4 : return (strength/(radius*radius*radius*12.5663706143591730))*evaluation3D_4thOrder(r2radius);
    case 6 : return (strength/(radius*radius*radius*12.5663706143591730))*evaluation3D_6thOrder(r2radius);
    }
    return 0.0;
	}


//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double,double,double )> getEvaluationPtr() const
	{
	std::function<double(double,double,double)> F = [this](double x,double y, double z) {return this->operator()(x,y,z);};
	return F;
	}


    // Internal evaluation routines for a normalized argument (r/radius)^2

    private: 

    // returns 4*pi*mollifier value

    double evaluation3D_2ndOrder(double r2) const
	{
    double u = 1.0-r2;
	switch(this->exponent-1)
    {
	case  0 : return   (15.0/2.0)*u;  break;
    case  1 : return   (105.0/8.0)*u*u; break;
    case  2 : return   (315.0/16.0)*u*u*u; break;
    case  3 : return   (3465.0/128.0)*u*u*u*u; break;
    case  4 : return   (9009.0/256.0)*u*u*u*u*u; break;
    case  5 : return   (45045.0/1024.0)*u*u*u*u*u*u; break;
    case  6 : return   (109395.0/2048.0)*u*u*u*u*u*u*u; break;
    case  7 : return   (2078505.0/32768.0)*u*u*u*u*u*u*u*u; break;
    case  8 : return   (4849845.0/65536.0)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return   (22309287.0/262144.0)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
	}

	// returns 4*pi*mollifier value

	double evaluation3D_4thOrder(double r2) const
	{
    double u = 1.0-r2;
	switch(this->exponent-1)
    {
	case  0 : return(-105.0/4.0)*(u- 9.0/4.0*u*u);       break;
    case  1 : return(-945.0/16.0)*((1.0-11.0/6.0*u)*u*u);   break;
    case  2 : return(-3465.0/32.0)*((1.0-13.0/8.0*u)*u*u*u);  break;
    case  3 : return(-45045.0/256.0)*((1.0-3.0/2.0*u)*u*u*u*u);  break;
    case  4 : return(-135135.0/512.0)*((1.0-17.0/12.0*u)*u*u*u*u*u); break;
    case  5 : return(-765765.0/2048.0)*((1.0-19.0/14.0*u)*u*u*u*u*u*u); break;
    case  6 : return(-2078505.0/4096.0)*((1.0-21.0/16.0*u)*u*u*u*u*u*u*u);  break;
    case  7 : return(-43648605.0/65536.0)*((1.0-23.0/18.0*u)*u*u*u*u*u*u*u*u);  break;
    case  8 : return(-111546435.0/131072.0)*((1.0-5.0/4.0*u)*u*u*u*u*u*u*u*u*u);  break;
    case  9 : return(-557732175.0/524288.0)*((1.0-27.0/22.0*u)*u*u*u*u*u*u*u*u*u*u);break;
    }
    return 0.0;
	}

    // returns 4*pi*mollifier value

    double evaluation3D_6thOrder(double r2) const
	{
    double u = 1.0-r2;
	switch(this->exponent-1)
    {
	case  0 : return(945.0/16.0)*(u +(-11.0/2.0+143.0/24.0*u)*u*u); break;
    case  1 : return(10395.0/64.0)*((1.0+(-13.0/3.0+65.0/16.0*u)*u)*u*u); break;
    case  2 : return(45045.0/128.0)*((1.0+(-15.0/4.0+51.0/16.0*u)*u)*u*u*u); break;
    case  3 : return(675675.0/1024.0)*((1.0+(-17.0/5.0+323.0/120.0*u)*u)*u*u*u*u);  break;
    case  4 : return(2297295.0/2048.0)*((1.0+(-19.0/6.0+19.0/8.0*u)*u)*u*u*u*u*u);  break;
    case  5 : return(14549535.0/8192.0)*((1.0+(-3.0+69.0/32.0*u)*u)*u*u*u*u*u*u);    break;
    case  6 : return(43648605.0/16384.0)*((1.0+(-23.0/8.0+575.0/288.0*u)*u)*u*u*u*u*u*u*u);break;
    case  7 : return(1003917915.0/262144.0)*((1.0+(-25.0/9.0+15.0/8.0*u)*u)*u*u*u*u*u*u*u*u);  break;
    case  8 : return(2788660875.0/524288.0)*((1.0+(-27.0/10.0+783.0/440.0*u)*u)*u*u*u*u*u*u*u*u*u);  break;
    case  9 : return(15058768725.0/2097152.0)*((1.0+(-29.0/11.0+899.0/528.0*u)*u)*u*u*u*u*u*u*u*u*u*u);  break;
    }
    return 0.0;
	}


    double        xPos;    // The position of the mollifier
    double        yPos;
    double        zPos;

    double      radius;    // The radius of the mollifier
    double    strength;    // The strength of the mollifier

    int           order;    // Order of the mollifier = 2, 4 or 6
    int        exponent;    // The exponent of the mollifier (determines differentiability)
};

#undef MOLLIFIER_DEFAULT_ORDER
#undef DEFAULT_DIFFERENTIABLITY


#endif /* SMOOTHPOLYMOLLIFIER_H_ */
