/*
 * HighOrderPolyMollifier1d.h
 *
 *  Created on: July 2, 2015
 *      Author: anderson
 *
 *  Updated July 2, 2015
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
#include <functional>

#ifndef _HighOrderPolyMollifier1d_
#define _HighOrderPolyMollifier1d_

#define _DEFAULT_ORDER_      4
#define _DEFAULT_DIFFERENTIABLITY_ 6

class HighOrderPolyMollifier1d
{
	public:

	HighOrderPolyMollifier1d()
	{
	initialize();
	}

	HighOrderPolyMollifier1d(const HighOrderPolyMollifier1d& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    radius     = S.radius;
    order      = S.order;
    exponent   = S.exponent;
	}

    HighOrderPolyMollifier1d(double xPos, double radius, double strength)
	{
    initialize(xPos,radius,strength);
	}

	void initialize()
	{
	xPos      = 0.0;
    radius    = 0.0;
    strength  = 0.0;

    order     = _DEFAULT_ORDER_;
    exponent  = _DEFAULT_DIFFERENTIABLITY_ + 1;
	}

	void initialize(double xPos, double radius, double strength)
	{
    this->strength  = strength;
    this->radius     = radius;
    this->xPos      = xPos;

    this->exponent  = _DEFAULT_DIFFERENTIABLITY_ + 1;
    this->order     = _DEFAULT_ORDER_;
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

	void setDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;
	}

	int getDifferentiablity()
	{
	return this->exponent-1;
	}

	/// Evaluation operator 1D

	double operator()(double x) const
	{
	double r2radius = ((x-xPos)*(x-xPos))/(radius*radius);
    if(r2radius >= 1.0) return 0.0;
    switch(order)
    {
    case 2 : return (strength/radius)*evaluation1D_2ndOrder(r2radius);
    case 4 : return (strength/radius)*evaluation1D_4thOrder(r2radius);
    case 6 : return (strength/radius)*evaluation1D_6thOrder(r2radius);
    }
    return 0.0;
	}

//  Returns a std::function that is bound to the evaluation operator of *this


	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return std::move(F);
	}



    //
    // Internal evaluation routines for a normalized argument (r/radius)^2
    //

    private: 

    double evaluation1D_2ndOrder(double r2) const
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

    double evaluation1D_4thOrder(double r2) const
    {
    double u = 1.0-r2;
    switch(this->exponent-1)
    {
	case  0 : return (-15.0/8.0)*(u-7.0/4.0*u*u); break;
    case  1 : return (-105.0/32.0)*(1.0-3.0/2.0*u)*u*u; break;
    case  2 : return (-315.0/64.0)*(1.0-11.0/8.0*u)*u*u*u; break;
    case  3 : return (-3465.0/512.0)*(1.0-13.0/10.0*u)*u*u*u*u; break;
    case  4 : return (-9009.0/1024.0)*(1.0-5.0/4.0*u)*u*u*u*u*u; break;
    case  5 : return (-45045.0/4096.0)*(1.0-17.0/14.0*u)*u*u*u*u*u*u; break;
    case  6 : return (-109395.0/8192.0)*(1.0-19.0/16.0*u)*u*u*u*u*u*u*u; break;
    case  7 : return (-2078505.0/131072.0)*(1.0-7.0/6.0*u)*u*u*u*u*u*u*u*u; break;
    case  8 : return (-4849845.0/262144.0)*(1.0-23.0/20.0*u)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return (-22309287.0/1048576.0)* (1.0-25.0/22.0*u)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
    }

    double evaluation1D_6thOrder(double r2) const
    {
    double u = 1.0-r2;
    switch(this->exponent-1)
    {
	case  0 : return (105.0/32.0)*(u+(-9.0/2.0+33.0/8.0*u)*u*u); break;
    case  1 : return (945.0/128.0)*(1.0+(-11.0/3.0+143.0/48.0*u)*u)*u*u; break;
    case  2 : return (3465.0/256.0)*(1.0+(-13.0/4.0+39.0/16.0*u)*u)*u*u*u; break;
    case  3 : return (45045.0/2048.0)*(1.0+(-3.0+17.0/8.0*u)*u)*u*u*u*u; break;
    case  4 : return (135135.0/4096.0)*(1.0+(-17.0/6.0+323.0/168.0*u)*u)*u*u*u*u*u; break;
    case  5 : return (765765.0/16384.0)*(1.0+(-19.0/7.0+57.0/32.0*u)*u)*u*u*u*u*u*u; break;
    case  6 : return (2078505.0/32768.0)*(1.0+(-21.0/8.0+161.0/96.0*u)*u)*u*u*u*u*u*u*u; break;
    case  7 : return (43648605.0/524288.0)* (1.0+(-23.0/9.0+115.0/72.0*u)*u)*u*u*u*u*u*u*u*u; break;
    case  8 : return (111546435.0/1048576.0)*(1.0+(-5.0/2.0+135.0/88.0*u)*u)*u*u*u*u*u*u*u*u*u; break;
    case  9 : return (557732175.0/4194304.0)*(1.0+(-27.0/11.0+261.0/176.0*u)*u)*u*u*u*u*u*u*u*u*u*u; break;
    }
    return 0.0;
    }

    double        xPos;    // The position of the mollifier
    double      radius;    // The radius of the mollifier
    double    strength;    // The strength of the mollifier

    int           order;    // Order of the mollifier = 2, 4 or 6
    int        exponent;    // The exponent of the mollifier (determines differentiability)
};

#undef _DEFAULT_ORDER_
#undef _DEFAULT_DIFFERENTIABLITY_


#endif /* SMOOTHPOLYMOLLIFIER_H_ */
