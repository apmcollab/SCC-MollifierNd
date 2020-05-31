/*
 * HighOrderPolyMollifier2d.h
 *
 *  Created on: Feb 16, 2014
 *      Author: anderson
 *
 *  Updated March 11, 2015
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

#ifndef HIGH_ORDER_POLY_MOLLIFIER_2D_
#define HIGH_ORDER_POLY_MOLLIFIER_2D_

#define MOLLIFIER_DEFAULT_ORDER      4
#define DEFAULT_DIFFERENTIABLITY 6

class HighOrderPolyMollifier2d
{
	public:

	HighOrderPolyMollifier2d()
	{
	initialize();
	}

	HighOrderPolyMollifier2d(const HighOrderPolyMollifier2d& S)
	{
	xPos       = S.xPos;
    yPos       = S.yPos;
    radius     = S.radius;
    strength   = S.strength;
    order      = S.order;
    exponent   = S.exponent;
	}



    HighOrderPolyMollifier2d(double xPos, double yPos, double radius, double strength)
	{
    initialize(xPos,yPos,radius, strength);
	}


	void initialize()
	{
    xPos      = 0.0;
    yPos      = 0.0;
    radius    = 0.0;
    strength  = 0.0;

    order     = MOLLIFIER_DEFAULT_ORDER;
    exponent  = DEFAULT_DIFFERENTIABLITY + 1;
	}

	void initialize(double xPos, double yPos, double radius, double strength)
	{
	this->xPos      = xPos;
    this->yPos      = yPos;
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
	{this->order = MOLLIFIER_DEFAULT_ORDER;}
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

    /// Evaluation operator 2D

	double operator()(double x, double y) const
	{
    double r2radius = ((x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos))/(radius*radius);
    if(r2radius >= 1.0) return 0.0;
    switch(order)
    {
    case 2 : return (strength/(radius*radius*6.28318530717958648))*evaluation2D_2ndOrder(r2radius);
    case 4 : return (strength/(radius*radius*6.28318530717958648))*evaluation2D_4thOrder(r2radius);
    case 6 : return (strength/(radius*radius*6.28318530717958648))*evaluation2D_6thOrder(r2radius);
    }
    return 0.0;
	}


//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double,double)> getEvaluationPtr() const
	{
	std::function<double(double,double)> F = [this](double x,double y) {return this->operator()(x,y);};
	return F;
	}



    //
    // Internal evaluation routines for a normalized argument (r/radius)^2
    //

    private: 


    // returns 2*pi*mollifier value

	double evaluation2D_2ndOrder(double r2) const
	{
    double u = 1.0-r2;

	switch(this->exponent-1)
    {
	case  0 : return   2*(exponent+1)*u;  break;
    case  1 : return   2*(exponent+1)*u*u; break;
    case  2 : return   2*(exponent+1)*u*u*u; break;
    case  3 : return   2*(exponent+1)*u*u*u*u; break;
    case  4 : return   2*(exponent+1)*u*u*u*u*u; break;
    case  5 : return   2*(exponent+1)*u*u*u*u*u*u; break;
    case  6 : return   2*(exponent+1)*u*u*u*u*u*u*u; break;
    case  7 : return   2*(exponent+1)*u*u*u*u*u*u*u*u; break;
    case  8 : return   2*(exponent+1)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return   2*(exponent+1)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
	}

    // returns 2*pi*mollifier value

	double evaluation2D_4thOrder(double r2) const
	{
    double u = 1.0-r2;
	switch(this->exponent-1)
    {
	case  0 : return (-12.0)*(u-2.0*u*u);           break;
    case  1 : return (-24.0)*((1.0-5.0/3.0*u)*u*u);   break;
    case  2 : return (-40.0)*((1.0-3.0/2.0*u)*u*u*u);  break;
    case  3 : return(-60.0)*((1.0-7.0/5.0*u)*u*u*u*u);   break;
    case  4 : return(-84.0)*((1.0-4.0/3.0*u)*u*u*u*u*u);    break;
    case  5 : return(-112.0)*((1.0-9.0/7.0*u)*u*u*u*u*u*u);   break;
    case  6 : return(-144.0)*((1.0-5.0/4.0*u)*u*u*u*u*u*u*u);  break;
    case  7 : return(-180.0)*((1.0-11.0/9.0*u)*u*u*u*u*u*u*u*u);  break;
    case  8 : return(-220.0)*((1.0-6.0/5.0*u)*u*u*u*u*u*u*u*u*u);    break;
    case  9 : return(-264.0)*((1.0-13.0/11.0*u)*u*u*u*u*u*u*u*u*u*u);  break;
    }
    return 0.0;
	}

    // returns 2*pi*mollifier value

    double evaluation2D_6thOrder(double r2) const
	{
    double u = 1.0-r2;
   	switch(this->exponent-1)
    {
	case  0 : return( 24.0)*(u+(-5.0+5.0*u)*u*u);  break;
    case  1 : return(60.0)*((1.0+(-4.0+7.0/2.0*u)*u)*u*u);  break;
    case  2 : return(120.0)*((1.0+(-7.0/2.0+14.0/5.0*u)*u)*u*u*u);  break;
    case  3 : return(210.0)*((1.0+(-16.0/5.0+12.0/5.0*u)*u)*u*u*u*u);  break;
    case  4 : return(336.0)*((1.0+(-3.0+15.0/7.0*u)*u)*u*u*u*u*u);       break;
    case  5 : return(504.0)*((1.0+(-20.0/7.0+55.0/28.0*u)*u)*u*u*u*u*u*u);  break;
    case  6 : return(720.0)*((1.0+(-11.0/4.0+11.0/6.0*u)*u)*u*u*u*u*u*u*u); break;
    case  7 : return(990.0)*((1.0+(-8.0/3.0+26.0/15.0*u)*u)*u*u*u*u*u*u*u*u);  break;
    case  8 : return(1320.0)*((1.0+(-13.0/5.0+91.0/55.0*u)*u)*u*u*u*u*u*u*u*u*u);  break;
    case  9 : return(1716.0)*((1.0+(-28.0/11.0+35.0/22.0*u)*u)*u*u*u*u*u*u*u*u*u*u);  break;
    }
    return 0.0;
	}

    double        xPos;    // The position of the mollifier
    double        yPos;
    double      radius;    // The radius of the mollifier
    double    strength;    // The strength of the mollifier

    int           order;    // Order of the mollifier = 2, 4 or 6
    int        exponent;    // The exponent of the mollifier (determines differentiability)
};

#undef MOLLIFIER_DEFAULT_ORDER
#undef DEFAULT_DIFFERENTIABLITY


#endif /* SMOOTHPOLYMOLLIFIER_H_ */
