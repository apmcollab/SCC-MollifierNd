/*
 * SymmetricPolyMollifier1d.h
 *
 *  Created on: Dec 26, 2023
 *      Author: anderson
 *
 * A class that determines a continuous function such that
 * the moments about xPos match up to the specified order
 * (2 through 8) those of a delta function located at xPos.
 *
 * Symmetric mollifiers of order p = 2..8 that are based on
 * the mollifer M(x) centered at x = xBar and of the form 
 *
 *
 * M(x) = (1/radius)*[ (1 - ((x-Bar)/radius)^2)^r * polyFactor((x-xBar)/radius) ]
 *
 * Where polyFactor is a polynomial defined by the relation
 *
 *                 k = p/2 -1
 * polyFactor(s) =SUM   a[k]*LegendrePoly(2*k,s)
 *                 k = 0
 *
 * in which LegendrePoly(2*k,s) is the 2*kth Legendre polynomial
 * with coefficients,a[k], determined so the requisite moments vanish.
 *
 *
 * This mollifer vanishes outside of [xBar - radius,xBar + radius]
 * and it's extension to [-oo, oo] by zero is r-1 times continuously
 * differentiable.
 *
 *
 * The coefficients a[k] are determined so that for p < order
 *
 *      -- radius 
 *     |                  | 1  if p = 0
 *     |  M(s)*s^p ds =  -|
 *     |                  | 0  if 0 < p < order
 *    -- -radius
 *
 * If the order specified is odd, then a mollfier of 
 * one higher order is implemented, since symmetry results 
 * in a higher order mollfier. 
 *
 */
 /*
#############################################################################
#
# Copyright 2023 - Chris Anderson
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
#include <vector>

#ifndef SYMMETRIC_POLY_MOLLIFIER_1D_
#define SYMMETRIC_POLY_MOLLIFIER_1D_

#define MOLLIFIER_DEFAULT_ORDER  4
#define DEFAULT_DIFFERENTIABLITY 6

#include "OrthoFunctionNd/SCC_PolyFun.h"

class SymmetricPolyMollifier1d
{

public:



	SymmetricPolyMollifier1d()
	{
	initialize();
	}

	SymmetricPolyMollifier1d(const SymmetricPolyMollifier1d& S)
	{
    initialize(S);
	}

    SymmetricPolyMollifier1d(double xPos, double radius, double strength)
	{
    initialize(xPos,radius,strength);
	}

	void initialize()
	{
	xPos      = 0.0;
    radius    = 0.0;
    strength  = 0.0;

    order     = MOLLIFIER_DEFAULT_ORDER;
    exponent  = DEFAULT_DIFFERENTIABLITY + 1;

    polyFactor.initialize();
    DpolyFactor.initialize();

    initializeCoefficients();
	}

	void initialize(const SymmetricPolyMollifier1d& S)
	{
    strength   = S.strength;
    xPos       = S.xPos;
    radius     = S.radius;
    order      = S.order;
    exponent   = S.exponent;

    initializeCoefficients();

    polyFactor.initialize(S.polyFactor);
    DpolyFactor.initialize(S.polyFactor);
	}

	void initialize(double xPos, double radius, double strength)
	{
		initialize(xPos, radius, strength, MOLLIFIER_DEFAULT_ORDER, DEFAULT_DIFFERENTIABLITY);
	}


	void initialize(double xPos, double radius, double strength,int order, int differentiability)
	{
    this->strength  = strength;
    this->radius    = radius;
    this->xPos      = xPos;

    this->exponent  = differentiability + 1;
    if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;

    this->order  = order;
    if(order%2 != 0)    {this->order = order+1;}
	if(this->order > 8) {this->order = 8;}

    polyFactor.initialize();
    DpolyFactor.initialize();

    initializeCoefficients();
	setPolyFactor();
	}

	void setRadius(double radius) {this->radius = radius;}
	double getRadius()            {return this->radius;  }

    void setStrength(double strength) {this->strength = strength;}
	double getStrength() const        {return this->strength;}

	void setParameters(int order, int differentiability)
	{
    // For odd order, bump up order by 1, since
    // it's no extra work due to symmetry.

    this->order = order;
	if(order%2 != 0){this->order = order+1;}

	if(this->order > 8) {this->order = 8;}

    this->exponent = differentiability + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;

    // Construct the polynomial multiplicative factor

    setPolyFactor();
	}


	void setOrder(int order)
	{
    // For odd order, bump up order by 1, since
    // it's no extra work due to symmetry.

	this->order = order;
	if(order%2 != 0){this->order = order+1;}

	if(this->order > 8) {this->order = 8;}

	setPolyFactor();
	}

	int getOrder() const
	{return this->order;}

	// 0 <= diffOrder <= 9

	void setDifferentiability(int diffOrder)
	{
	this->exponent = diffOrder + 1;
	if(exponent > 10) this->exponent = 10;
    if(exponent < 1 ) this->exponent = 1;
    setPolyFactor();
	}

	int getDifferentiablity()
	{
	return this->exponent-1;
	}

	void setPolyFactor()
	{
	polyFactor.initialize(order - 2);
    long k = exponent-1;
	switch(order)
    {
	case  2 : polyFactor[0] = a2_0[k]; break;

    case  4 : polyFactor[0] = a4_0[k];
              polyFactor[2] = a4_2[k]; break;

    case  6 : polyFactor[0] = a6_0[k];
              polyFactor[2] = a6_2[k];
              polyFactor[4] = a6_4[k]; break;

    case  8 : polyFactor[0] = a8_0[k];
              polyFactor[2] = a8_2[k];
              polyFactor[4] = a8_4[k];
              polyFactor[6] = a8_6[k]; break;
    }

	/*
	//
	// Code based on expansion in Legendre polynomials
	//
	SCC::OrthoPoly orthogonalPoly(SCC::OrthoPoly::Legendre);
	long k = exponent-1;
	switch(order)
    {
	case  2 : polyFactor[0] = a2_0[k]; break;

    case  4 : polyFactor[0] = a4_0[k];
              polyFactor    = polyFactor + a4_2[k]*orthogonalPoly.getNthOrthoPolynomial(2); break;

    case  6 : polyFactor[0] = a6_0[k];
              polyFactor    = polyFactor + a6_2[k]*orthogonalPoly.getNthOrthoPolynomial(2);
              polyFactor    = polyFactor + a6_4[k]*orthogonalPoly.getNthOrthoPolynomial(4); break;

    case  8 : polyFactor[0] = a8_0[k];
              polyFactor    = polyFactor + a8_2[k]*orthogonalPoly.getNthOrthoPolynomial(2);
              polyFactor    = polyFactor + a8_4[k]*orthogonalPoly.getNthOrthoPolynomial(4);
              polyFactor    = polyFactor + a8_6[k]*orthogonalPoly.getNthOrthoPolynomial(6); break;

    }
    */

	DpolyFactor.initialize();
	DpolyFactor = polyFactor.differentiate();
	}

	double operator()(double x) const
	{
	double r2 = ((x-xPos)*(x-xPos))/(radius*radius);
    if(r2 >= 1.0) {return 0.0;}

    double u       = 1.0-r2;
    double pFactor = pow(u,exponent);

    double r       = (x-xPos)/radius;


    return (strength/radius)*pFactor*polyFactor(r);

    /*
    // Code based on expansion in Legendre polynomials
    long k         = exponent -1;
    switch(order)
    {
	case  2 : return  (strength/radius)*pFactor*(a2_0[k]); break;
    case  4 : return  (strength/radius)*pFactor*(a4_0[k] + a4_2[k]*0.5*(3.0*r2 - 1.0)); break;
    case  6 : return  (strength/radius)*pFactor*(a6_0[k] + a6_2[k]*0.5*(3.0*r2 - 1.0) + a6_4[k]*.125*(35.0*r2*r2 - 30.0*r2 + 3.0)); break;
    case  8 : return  (strength/radius)*pFactor*(a8_0[k] + a8_2[k]*0.5*(3.0*r2 - 1.0) + a8_4[k]*.125*(35.0*r2*r2 - 30.0*r2 + 3.0)
    		                                   + a8_6[k]*.0625*(231.0*r2*r2*r2 - 315.0*r2*r2 + 105.0*r2 - 5.0)); break;
    }
    return 0.0;
    */
	}

	//  Returns a std::function that is bound to the evaluation operator of *this

	std::function<double(double)> getEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->operator()(x);};
	return F;
	}

	double derivative(double x) const
    {
	double r2 = ((x-xPos)*(x-xPos))/(radius*radius);

    if(r2 >= 1.0) {return 0.0;}

    double u         = 1.0-r2;
    double pFactorM1 = pow(u,exponent-1);
    double pFactor   = pFactorM1*u;


    double r2Factor  = (2.0*(x-xPos))/(radius*radius);
    double rFactor   = 1.0/radius;
    double r         = (x-xPos)/radius;


    return (strength/radius)*(pFactorM1*exponent*(-1.0)*polyFactor(r)*r2Factor  + pFactor*DpolyFactor(r)*rFactor);

    /*
    // Code based on expansion in Legendre polynomials
    long k           = exponent -1;
    switch(order)
    {
	case  2 : return   (strength/radius)*pFactorM1*exponent*a2_0[k]*(-1.0)*r2Factor; break;

    case  4 : return  (strength/radius)*r2Factor*((-1.0)*pFactorM1*exponent*(a4_0[k]
																		  + a4_2[k]*0.5*(3.0*r2 - 1.0))
    		                                                     + pFactor*(a4_2[k]*1.5)); break;

    case  6 : return  (strength/radius)*r2Factor*((-1.0)*pFactorM1*exponent*(a6_0[k]
																		  + a6_2[k]*0.5*(3.0*r2 - 1.0)
																		  + a6_4[k]*.125*(35.0*r2*r2 - 30.0*r2 + 3.0))
    		                                                     + pFactor*(a6_2[k]*1.5
    		                                                    		  + a6_4[k]*.125*(70.0*r2 - 30.0))); break;

    case  8 : return  (strength/radius)*r2Factor*((-1.0)*pFactorM1*exponent*(a8_0[k] + a8_2[k]*0.5*(3.0*r2 - 1.0)
    		                                                              + a8_4[k]*.125*(35.0*r2*r2 - 30.0*r2 + 3.0)
    	                                                                  + a8_6[k]*.0625*(231.0*r2*r2*r2 - 315.0*r2*r2 + 105.0*r2 - 5.0))
			                                                     + pFactor*(a8_2[k]*1.5
			                                                    		  + a8_4[k]*.125*(70.0*r2  - 30.0)
			                                                    		  + a8_6[k]*.0625*(693.0*r2*r2 - 630.0*r2 + 105.0))); break;
    }
    return  0.0;
    */
    }


	//  Returns a std::function that is bound to the derivative evaluation operator of *this

	std::function<double(double)> getDerivativeEvaluationPtr() const
	{
	std::function<double(double)> F = [this](double x) {return this->derivative(x);};
	return F;
	}

	void initializeCoefficients()
	{

	a2_0.resize(10,0.0);

	a4_0.resize(10,0.0);
	a4_2.resize(10,0.0);

	a6_0.resize(10,0.0);
	a6_2.resize(10,0.0);
	a6_4.resize(10,0.0);

    a8_0.resize(10,0.0);
	a8_2.resize(10,0.0);
	a8_4.resize(10,0.0);
	a8_6.resize(10,0.0);

/* 2th order coefficients */

	a2_0[0] = 3.0/4.0;
	a2_0[1] = 15.0/16.0;
	a2_0[2] = 35.0/32.0;
	a2_0[3] = 315.0/256.0;
	a2_0[4] = 693.0/512.0;
	a2_0[5] = 3003.0/2048.0;
	a2_0[6] = 6435.0/4096.0;
	a2_0[7] = 109395.0/65536.0;
	a2_0[8] = 230945.0/131072.0;
	a2_0[9] = 969969.0/524288.0;

	// 4th order polynomial  coefficients
	
	a4_0[0] = 45.0/32.0;
	a4_0[1] = 105.0/64.0;
	a4_0[2] = 945.0/512.0;
	a4_0[3] = 2079.0/1024.0;
	a4_0[4] = 9009.0/4096.0;
	a4_0[5] = 19305.0/8192.0;
	a4_0[6] = 328185.0/131072.0;
	a4_0[7] = 692835.0/262144.0;
	a4_0[8] = 2909907.0/1048576.0;
	a4_0[9] = 6084351.0/2097152.0;

    a4_2[0] = -105.0/32.0;
    a4_2[1] = -315.0/64.0;
    a4_2[2] = -3465.0/512.0;
    a4_2[3] = -9009.0/1024.0;
    a4_2[4] = -45045.0/4096.0;
    a4_2[5] = -109395.0/8192.0;
    a4_2[6] = -2078505.0/131072.0;
    a4_2[7] = -4849845.0/262144.0;
    a4_2[8] = -22309287.0/1048576.0;
    a4_2[9] = -50702925.0/2097152.0;

    // 6th order polynomial coefficients

    a6_0[0] = 525.0/256.0;
    a6_0[1] = 4725.0/2048.0;
    a6_0[2] = 10395.0/4096.0;
    a6_0[3] = 45045.0/16384.0;
    a6_0[4] = 96525.0/32768.0;
    a6_0[5] = 1640925.0/524288.0;
    a6_0[6] = 3464175.0/1048576.0;
    a6_0[7] = 14549535.0/4194304.0;
    a6_0[8] = 30421755.0/8388608.0;
    a6_0[9] = 253514625.0/67108864.0;

    a6_2[0] = -1575.0/128.0;
    a6_2[1] = -17325.0/1024.0;
    a6_2[2] = -45045.0/2048.0;
    a6_2[3] = -225225.0/8192.0;
    a6_2[4] = -546975.0/16384.0;
    a6_2[5] = -10392525.0/262144.0;
    a6_2[6] = -24249225.0/524288.0;
    a6_2[7] = -111546435.0/2097152.0;
    a6_2[8] = -253514625.0/4194304.0;
    a6_2[9] = -2281631625.0/33554432.0;

    a6_4[0] = 3465.0/256.0;
    a6_4[1] = 45045.0/2048.0;
    a6_4[2] = 135135.0/4096.0;
    a6_4[3] = 765765.0/16384.0;
    a6_4[4] = 2078505.0/32768.0;
    a6_4[5] = 43648605.0/524288.0;
    a6_4[6] = 111546435.0/1048576.0;
    a6_4[7] = 557732175.0/4194304.0;
    a6_4[8] = 1368978975.0/8388608.0;
    a6_4[9] = 13233463425.0/67108864.0;


//  8th order polynomial coefficients
	
    a8_0[0] = 11025.0/4096.0;
    a8_0[1] = 24255.0/8192.0;
    a8_0[2] = 105105.0/32768.0;
    a8_0[3] = 225225.0/65536.0;
    a8_0[4] = 3828825.0/1048576.0;
    a8_0[5] = 8083075.0/2097152.0;
    a8_0[6] = 33948915.0/8388608.0;
    a8_0[7] = 70984095.0/16777216.0;
    a8_0[8] = 591534125.0/134217728.0;
    a8_0[9] = 1228570875.0/268435456.0;

    a8_2[0] = -121275.0/4096.0;
    a8_2[1] = -315315.0/8192.0;
    a8_2[2] = -1576575.0/32768.0;
    a8_2[3] = -3828825.0/65536.0;
    a8_2[4] = -72747675.0/1048576.0;
    a8_2[5] = -169744575.0/2097152.0;
    a8_2[6] = -780825045.0/8388608.0;
    a8_2[7] = -1774602375.0/16777216.0;
    a8_2[8] = -15971421375.0/134217728.0;
    a8_2[9] = -35628555375.0/268435456.0;

    a8_4[0] = 315315.0/4096.0;
    a8_4[1] = 945945.0/8192.0;
    a8_4[2] = 5360355.0/32768.0;
    a8_4[3] = 14549535.0/65536.0;
    a8_4[4] = 305540235.0/1048576.0;
    a8_4[5] = 780825045.0/2097152.0;
    a8_4[6] = 3904125225.0/8388608.0;
    a8_4[7] = 9582852825.0/16777216.0;
    a8_4[8] = 92634243975.0/134217728.0;
    a8_4[9] = 220897043325.0/268435456.0;

    a8_6[0] = -225225.0/4096.0;
    a8_6[1] = -765765.0/8192.0;
    a8_6[2] = -4849845.0/32768.0;
    a8_6[3] = -14549535.0/65536.0;
    a8_6[4] = -334639305.0/1048576.0;
    a8_6[5] = -929553625.0/2097152.0;
    a8_6[6] = -5019589575.0/8388608.0;
    a8_6[7] = -13233463425.0/16777216.0;
    a8_6[8] = -136745788725.0/134217728.0;
    a8_6[9] = -347123925225.0/268435456.0;
/*
// Coefficients for code based on expansion in Legendre polynomials
4th order Legendre coefficients  

      a4_0[0] = 5.0/16.0;
      a4_0[1] = 0.0;
      a4_0[2] = -105.0/256.0;
      a4_0[3] = -231.0/256.0;
      a4_0[4] = -3003.0/2048.0;
      a4_0[5] = -2145.0/1024.0;
      a4_0[6] = -182325.0/65536.0;
      a4_0[7] = -230945.0/65536.0;
      a4_0[8] = -2263261.0/524288.0;
      a4_0[9] = -676039.0/131072.0;

      a4_2[0] = -35.0/16.0;
      a4_2[1] = -105.0/32.0;
      a4_2[2] = -1155.0/256.0;
      a4_2[3] = -3003.0/512.0;
      a4_2[4] = -15015.0/2048.0;
      a4_2[5] = -36465.0/4096.0;
      a4_2[6] = -692835.0/65536.0;
      a4_2[7] = -1616615.0/131072.0;
      a4_2[8] = -7436429.0/524288.0;
      a4_2[9] = -16900975.0/1048576.0;

6th order Legendre coefficients 

      a6_0[0] = 21.0/32.0;
      a6_0[1] = 273.0/256.0;
      a6_0[2] = 231.0/128.0;
      a6_0[3] = 3003.0/1024.0;
      a6_0[4] = 18447.0/4096.0;
      a6_0[5] = 430287.0/65536.0;
      a6_0[6] = 600457.0/65536.0;
      a6_0[7] = 1616615.0/131072.0;
      a6_0[8] = 16900975.0/1048576.0;
      a6_0[9] = 172389945.0/8388608.0;

      a6_2[0] = -15.0/32.0;
      a6_2[1] = 165.0/128.0;
      a6_2[2] = 2145.0/512.0;
      a6_2[3] = 2145.0/256.0;
      a6_2[4] = 401115.0/28672.0;
      a6_2[5] = 692835.0/32768.0;
      a6_2[6] = 3926065.0/131072.0;
      a6_2[7] = 5311735.0/131072.0;
      a6_2[8] = 55531775.0/1048576.0;
      a6_2[9] = 282487725.0/4194304.0;

      a6_4[0] = 99.0/32.0;
      a6_4[1] = 1287.0/256.0;
      a6_4[2] = 3861.0/512.0;
      a6_4[3] = 21879.0/2048.0;
      a6_4[4] = 415701.0/28672.0;
      a6_4[5] = 1247103.0/65536.0;
      a6_4[6] = 3187041.0/131072.0;
      a6_4[7] = 15935205.0/524288.0;
      a6_4[8] = 39113685.0/1048576.0;
      a6_4[9] = 378098955.0/8388608.0;


8th order Legendre coefficients 

      a8_0[0] = 93.0/256.0;
      a8_0[1] = -33.0/256.0;
      a8_0[2] = -1287.0/1024.0;
      a8_0[3] = -429.0/128.0;
      a8_0[4] = -444873.0/65536.0;
      a8_0[5] = -785213.0/65536.0;
      a8_0[6] = -2540395.0/131072.0;
      a8_0[7] = -482885.0/16384.0;
      a8_0[8] = -358783555.0/8388608.0;
      a8_0[9] = -501791805.0/8388608.0;

      a8_2[0] = -495.0/256.0;
      a8_2[1] = -2145.0/512.0;
      a8_2[2] = -9295.0/1024.0;
      a8_2[3] = -36465.0/2048.0;
      a8_2[4] = -2078505.0/65536.0;
      a8_2[5] = -20554105.0/393216.0;
      a8_2[6] = -5311735.0/65536.0;
      a8_2[7] = -31387525.0/262144.0;
      a8_2[8] = -1426925175.0/8388608.0;
      a8_2[9] = -3926412225.0/16777216.0;

      a8_4[0] = 117.0/256.0;
      a8_4[1] = -351.0/128.0;
      a8_4[2] = -17901.0/2048.0;
      a8_4[3] = -37791.0/2048.0;
      a8_4[4] = -2154087.0/65536.0;
      a8_4[5] = -869193.0/16384.0;
      a8_4[6] = -42010995.0/524288.0;
      a8_4[7] = -664932645.0/5767168.0;
      a8_4[8] = -14745859245.0/92274688.0;
      a8_4[9] = -901620585.0/4194304.0;

      a8_6[0] = -975.0/256.0;
      a8_6[1] = -3315.0/512.0;
      a8_6[2] = -20995.0/2048.0;
      a8_6[3] = -62985.0/4096.0;
      a8_6[4] = -1448655.0/65536.0;
      a8_6[5] = -12072125.0/393216.0;
      a8_6[6] = -21729825.0/524288.0;
      a8_6[7] = -630164925.0/11534336.0;
      a8_6[8] = -6511704225.0/92274688.0;
      a8_6[9] = -1502700975.0/16777216.0;
      */
	}


    std::vector<double> a2_0; // second order coefficients

	std::vector<double> a4_0; // fourth order coefficients
	std::vector<double> a4_2;

	std::vector<double> a6_0; // sixth order coefficients
	std::vector<double> a6_2;
	std::vector<double> a6_4;

    std::vector<double> a8_0; // eight order coefficients
	std::vector<double> a8_2;
	std::vector<double> a8_4;
	std::vector<double> a8_6;

	double        xPos;    // The position of the mollifier
    double      radius;    // The radius of the mollifier
    double    strength;    // The strength of the mollifier

    int           order;    // Order of the mollifier = 2-8 (odd order specification results in next higher order)
    int        exponent;    // The exponent of the mollifier (determines differentiability)

    SCC::PolyFun  polyFactor;
    SCC::PolyFun DpolyFactor;
};

#undef MOLLIFIER_DEFAULT_ORDER
#undef DEFAULT_DIFFERENTIABLITY

#endif /* SYMMETRIC_POLY_MOLLIFIER_1D_ */
