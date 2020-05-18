/*
 * HighOrderPolyPotential3d.h
 *
 *  Created on: July 3, 2015
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

#ifndef _HighOrderPolyPotential3d_
#define _HighOrderPolyPotential3d_

#define _DEFAULT_ORDER_      4
#define _DEFAULT_SOURCE_DIFFERENTIABLITY_ 6

class HighOrderPolyPotential3d
{
	public:

	HighOrderPolyPotential3d()
	{
	initialize();
	}

	HighOrderPolyPotential3d(const HighOrderPolyPotential3d& S)
	{
    this->xPos       = S.xPos;
    this->yPos       = S.yPos;
    this->zPos       = S.zPos;
    this->radius     = S.radius;
    this->laplaceCoeff       = S.laplaceCoeff;
    this->strength   = S.strength;
    this->order      = S.order;
    this->exponent   = S.exponent;
	}



    HighOrderPolyPotential3d(double xPos, double yPos,double zPos, double radius, double strength,double laplaceCoeff = 1.0)
	{
    initialize(xPos,yPos,zPos, radius, strength, laplaceCoeff);
	}

	void initialize()
	{
    xPos      = 0.0;
    yPos      = 0.0;
    zPos      = 0.0;
    radius    = 0.0;
    laplaceCoeff      = 1.0;
    strength  = 0.0;

    order     = _DEFAULT_ORDER_;
    exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
	}



	void initialize( double xPos, double yPos, double zPos, double radius, double strength,double laplaceCoeff = 1.0)
	{

    this->xPos      = xPos;
    this->yPos      = yPos;
    this->zPos      = zPos;

    this->laplaceCoeff      = laplaceCoeff;
    this->radius    = radius;
    this->strength  = strength;

    this->exponent  = _DEFAULT_SOURCE_DIFFERENTIABLITY_ + 1;
    this->order     = _DEFAULT_ORDER_;
	}

	void setRadius(double radius) {this->radius = radius;}
	double getRadius()            {return this->radius;}

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



	double operator()(double x, double y, double z) const
	{
	double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos);
	double r2radius = r2/(radius*radius);
	if(r2radius >= 1.0)  {return -(strength/laplaceCoeff)*0.79577471545948e-01*(1.0/(std::sqrt(r2)));} // strength*(-1.0/(4*pi*r))
    switch(order)
    {
    case 2 : return ((evaluation3D_2ndOrder(r2radius) - 1.0)*.79577471545947667883e-01*(strength/laplaceCoeff))/radius;
    case 4 : return ((evaluation3D_4thOrder(r2radius) - 1.0)*.79577471545947667883e-01*(strength/laplaceCoeff))/radius;
    case 6 : return ((evaluation3D_6thOrder(r2radius) - 1.0)*.79577471545947667883e-01*(strength/laplaceCoeff))/radius;
    }
    return 0.0;
	}


	//  Returns a std::function that is bound to the evaluation operator of *this

 	std::function<double(double,double,double )> getEvaluationPtr() const
	{
	std::function<double(double,double,double)> F = [this](double x,double y, double z) {return this->operator()(x,y,z);};
	return F;
	}

    double radialDerivative3D(double r) const
    {
	double  dVal;
	double r2radius = (r*r)/(radius*radius);

	if(r2radius <= 1.0)
    {
	switch(order)
	{
    	case 2 : dVal = evaluate2xDq3D_2ndOrder(r2radius)*(r/radius); break;
    	case 4 : dVal = evaluate2xDq3D_4thOrder(r2radius)*(r/radius); break;
    	case 6 : dVal = evaluate2xDq3D_6thOrder(r2radius)*(r/radius); break;
	}
	return (dVal*(strength/laplaceCoeff)*.79577471545947667883e-01)/(radius*radius);
    }
    else
    {
    return (strength/laplaceCoeff)*.79577471545947667883e-01*(1.0/(r*r));
    }
    return 0.0;
    }

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =   z derivative


	void derivatives(double x, double y, double z, std::vector <double>& derivativeList, int maxOrder = 1) const
   	{
   	switch (maxOrder)
   	{
   	case 0   :  derivativeList.resize(1); break;
   	case 1   :  derivativeList.resize(4); break;
   	default  :  derivativeList.resize(4); break;
   	}

   	double xRadius  = (x-xPos)/radius;
   	double yRadius  = (y-yPos)/radius;
    double zRadius  = (z-zPos)/radius;
    double r2 = (x-xPos)*(x-xPos)  + (y-yPos)*(y-yPos) + (z-zPos)*(z-zPos) ;
	double r2radius = r2/(radius*radius);
	double r;
	double r3;

    int dCount = derivativeList.size();

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

    return;
   	}

    private :

    // returns 4*pi*mollifier value

    double evaluation3D_2ndOrder(double r2) const
	{
	switch(this->exponent-1)
	{
	case  0 : return    (15.0/2.0)*((1.0/6.0-r2/20.0)*r2 -  7.0/60.0); break;
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

    // Returns 4*Pi*2*q'(r^2) where u(r) = q(r^2)

    double evaluate2xDq3D_2ndOrder(double r2) const
	{
	switch(this->exponent-1)
	{
	case  0 : return    (15.0/2.0)*(1.0/3.0-r2/5.0); break;
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

	// returns 4*pi* value

	double evaluation3D_4thOrder(double r2) const
	{
	switch(this->exponent-1)
    {
	case  0 : return(-105.0/4.0)*(73.0/840.0+(-5.0/24.0+(7.0/40.0-3.0/56.0*r2)*r2)*r2);       break;
    case  1 : return(-945.0/16.0)*(689.0/15120.0+(-5.0/36.0+(7.0/40.0+(-3.0/28.0+11.0/432.0*r2)*r2)*r2)*r2);   break;
    case  2 : return(-3465.0/32.0)*(1567.0/55440.0+(-5.0/48.0+(7.0/40.0+(-9.0/56.0+(11.0/144.0-13.0/880.0*r2)*r2)*r2)*r2)*r2);  break;
    case  3 : return(-45045.0/256.0)*(6961.0/360360.0+(-1.0/12.0+(7.0/40.0+(-3.0/14.0+(11.0/72.0+(-13.0/220.0+r2/104)*r2)*r2)*r2)*r2)*r2);  break;
    case  4 : return(-135135.0/512.0)*(15209.0/1081080.0+(-5.0/72.0+(7.0/40.0+(-15.0/56.0+(55.0/216.0+(-13.0/88.0+(5.0/104.0-17.0/2520.0*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
    case  5 : return(-765765.0/2048.0)*(262649.0/24504480.0+(-5.0/84.0+(7.0/40.0+(-9.0/28.0+(55.0/144.0+(-13.0/44.0+(15.0/104.0+(-17.0/420.0+19.0/3808.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2); break;
    case  6 : return(-2078505.0/4096.0)*(561763.0/66512160.0+(-5.0/96.0+(7.0/40.0+(-3.0/8.0+(77.0/144.0+(-91.0/176.0+(35.0/104.0+(-17.0/120.0+(19.0/544.0-7.0/1824.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  7 : return(-43648605.0/65536.0)*(2385619.0/349188840.0+(-5.0/108.0+(7.0/40.0+(-3.0/7.0+(77.0/108.0+(-91.0/110.0+(35.0/52.0+(-17.0/45.0+(19.0/136.0+(-7.0/228.0+23.0/7560.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  8 : return(-111546435.0/131072.0)*(1007155.0/178474296.0+(-1.0/24.0+(7.0/40.0+(-27.0/56.0+(11.0/12.0+(-273.0/220.0+(63.0/52.0+(-17.0/20.0+(57.0/136.0+(-21.0/152.0+(23.0/840.0-5.0/2024.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  9 : return(-557732175.0/524288.0)*(42314317.0/8923714800.0+(-5.0/132.0+(7.0/40.0+(-15.0/28.0+(55.0/48.0+(-39.0/22.0+(105.0/52.0+(-17.0/10.0+(285.0/272.0+(-35.0/76.0+(23.0/168.0+(-25.0/1012.0+9.0/4400.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    }
    return 0.0;
	}

	// Returns 4*Pi*2*q'(r^2) where u(r) = q(r^2)

	double evaluate2xDq3D_4thOrder(double r2) const
	{
	switch(this->exponent-1)
    {
	case  0 : return 175.0/16.0+(-147.0/8.0+135.0/16.0*r2)*r2; break;
    case  1 : return 525.0/32.0+(-1323.0/32.0+(1215.0/32.0-385.0/32.0*r2)*r2)*r2;  break;
    case  2 : return 5775.0/256.0+(-4851.0/64.0+(13365.0/128.0+(-4235.0/64.0+4095.0/256.0*r2)*r2)*r2)*r2;  break;
    case  3 : return 15015.0/512.0+(-63063.0/512.0+(57915.0/256.0+(-55055.0/256.0+(53235.0/512.0-10395.0/512.0*r2)*r2)*r2)*r2)*r2;  break;
    case  4 : return 75075.0/2048.0+(-189189.0/1024.0+(868725.0/2048.0+(-275275.0/512.0+(798525.0/2048.0+(-155925.0/1024.0+51051.0/2048.0*r2)*r2)*r2)*r2)*r2)*r2; break;
    case  5 : return 182325.0/4096.0+(-1072071.0/4096.0+(2953665.0/4096.0+(-4679675.0/4096.0+(4524975.0/4096.0+(-2650725.0/4096.0+(867867.0/4096.0-122265.0/4096.0*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
    case  6 : return 3464175.0/65536.0+(-2909907.0/8192.0+(18706545.0/16384.0+(-17782765.0/8192.0+(85974525.0/32768.0+(-16787925.0/8192.0+(16489473.0/16384.0+(-2323035.0/8192.0+2297295.0/65536.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
    case  7 : return 8083075.0/131072.0+(-61108047.0/131072.0+(56119635.0/32768.0+(-124479355.0/32768.0+(361093005.0/65536.0+(-352546425.0/65536.0+(115426311.0/32768.0+(-48783735.0/32768.0+(48243195.0/131072.0-5311735.0/131072.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2;  break;
    case  8 : return 37182145.0/524288.0+(-156165009.0/262144.0+(1290751605.0/524288.0+(-409003595.0/65536.0+(2768379705.0/262144.0+(-1621713555.0/131072.0+(2654805153.0/262144.0+(-374008635.0/65536.0+(1109593485.0/524288.0+(-122169905.0/262144.0+24249225.0/524288.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
    case  9 : return 84504875.0/1048576.0+(-780825045.0/1048576.0+(3585421125.0/1048576.0+(-10225089875.0/1048576.0+(9887070375.0/524288.0+(-13514279625.0/524288.0+(13274025765.0/524288.0+(-9350215875.0/524288.0+(9246612375.0/1048576.0+(-3054247625.0/1048576.0+(606230625.0/1048576.0-54759159.0/1048576.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
    }
    return 0.0;
	}

    // returns 4*pi*mollifier value

    double evaluation3D_6thOrder(double r2) const
	{
	switch(this->exponent-1)
    {
	case  0 : return(945.0/16.0)*(-3701.0/60480.0+(35.0/144.0+(-63.0/160.0+(33.0/112.0-143.0/1728.0*r2)*r2)*r2)*r2); break;
    case  1 : return(10395.0/64.0)*(-8347.0/332640.0+(35.0/288.0+(-21.0/80.0+(33.0/112.0+(-143.0/864.0+13.0/352.0*r2)*r2)*r2)*r2)*r2); break;
    case  2 : return(45045.0/128.0)*(-36853.0/2882880.0+(7.0/96.0+(-63.0/320.0+(33.0/112.0+(-143.0/576.0+(39.0/352.0-17.0/832.0*r2)*r2)*r2)*r2)*r2)*r2); break;
    case  3 : return(675675.0/1024.0)*(-80141.0/10810800.0+(7.0/144.0+(-63.0/400.0+(33.0/112.0+(-143.0/432.0+(39.0/176.0+(-17.0/208.0+323.0/25200.0*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  4 : return(2297295.0/2048.0)*(-1378781.0/294053760.0+(5.0/144.0+(-21.0/160.0+(33.0/112.0+(-715.0/1728.0+(65.0/176.0+(-85.0/416.0+(323.0/5040.0-19.0/2176.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  5 : return(14549535.0/8192.0)*(-2939887.0/931170240.0+(5.0/192.0+(-9.0/80.0+(33.0/112.0+(-143.0/288.0+(195.0/352.0+(-85.0/208.0+(323.0/1680.0+(-57.0/1088.0+23.0/3648.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);    break;
    case  6 : return(43648605.0/16384.0)*(-12452383.0/5587021440.0+(35.0/1728.0+(-63.0/640.0+(33.0/112.0+(-1001.0/1728.0+(273.0/352.0+(-595.0/832.0+(323.0/720.0+(-399.0/2176.0+(161.0/3648.0-115.0/24192.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);break;
    case  7 : return(1003917915.0/262144.0)*(-26227451.0/16062686640.0+(7.0/432.0+(-7.0/80.0+(33.0/112.0+(-143.0/216.0+(91.0/88.0+(-119.0/104.0+(323.0/360.0+(-133.0/272.0+(161.0/912.0+(-115.0/3024.0+15.0/4048.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  8 : return(2788660875.0/524288.0)*(-219960193.0/178474296000.0+(7.0/528.0+(-63.0/800.0+(33.0/112.0+(-143.0/192.0+(117.0/88.0+(-357.0/208.0+(323.0/200.0+(-1197.0/1088.0+(161.0/304.0+(-115.0/672.0+(135.0/4048.0-261.0/88000.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    case  9 : return(15058768725.0/2097152.0)*(-459421511.0/481880599200.0+(35.0/3168.0+(-63.0/880.0+(33.0/112.0+(-715.0/864.0+(585.0/352.0+(-255.0/104.0+(323.0/120.0+(-1197.0/544.0+(805.0/608.0+(-575.0/1008.0+(675.0/4048.0+(-261.0/8800.0+899.0/370656.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2);  break;
    }
    return 0.0;
	}

    // Returns 4*Pi*2*q'(r^2) where u(r) = q(r^2)

    double evaluate2xDq3D_6thOrder(double r2) const
	{
	switch(this->exponent-1)
	{
    case  0 : return    3675.0/128.0+(-11907.0/128.0+(13365.0/128.0-5005.0/128.0*r2)*r2)*r2; break;
    case  1 : return    40425.0/1024.0+(-43659.0/256.0+(147015.0/512.0+(-55055.0/256.0+61425.0/1024.0*r2)*r2)*r2)*r2; break;
    case  2 : return    105105.0/2048.0+(-567567.0/2048.0+(637065.0/1024.0+(-715715.0/1024.0+(798525.0/2048.0-176715.0/2048.0*r2)*r2)*r2)*r2)*r2; break;
    case  3 : return    525525.0/8192.0+(-1702701.0/4096.0+(9555975.0/8192.0+(-3578575.0/2048.0+(11977875.0/8192.0+(-2650725.0/4096.0+969969.0/8192.0*r2)*r2)*r2)*r2)*r2)*r2;  break;
    case  4 : return    1276275.0/16384.0+(-9648639.0/16384.0+(32490315.0/16384.0+(-60835775.0/16384.0+(67874625.0/16384.0+(-45062325.0/16384.0+(16489473.0/16384.0-2567565.0/16384.0*r2)*r2)*r2)*r2)*r2)*r2)*r2;  break;
    case  5 : return    24249225.0/262144.0+(-26189163.0/32768.0+(205771995.0/65536.0+(-231175945.0/32768.0+(1289617875.0/131072.0+(-285394725.0/32768.0+(313299987.0/65536.0+(-48783735.0/32768.0+52837785.0/262144.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2;    break;
    case  6 : return    56581525.0/524288.0+(-549972423.0/524288.0+(617315985.0/131072.0+(-1618231615.0/131072.0+(5416395075.0/262144.0+(-5993289225.0/262144.0+(2193099909.0/131072.0+(-1024458435.0/131072.0+(1109593485.0/524288.0-132793375.0/524288.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  7 : return    260275015.0/2097152.0+(-1405485081.0/1048576.0+(14198267655.0/2097152.0+(-5317046735.0/262144.0+(41525695575.0/1048576.0+(-27569130435.0/524288.0+(50441297907.0/1048576.0+(-7854181335.0/262144.0+(25520650155.0/2097152.0+(-3054247625.0/1048576.0+654729075.0/2097152.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  8 : return    591534125.0/4194304.0+(-7027425405.0/4194304.0+(39439632375.0/4194304.0+(-132926168375.0/4194304.0+(148306055625.0/2097152.0+(-229742753625.0/2097152.0+(252206489535.0/2097152.0+(-196354533375.0/2097152.0+(212672084625.0/4194304.0+(-76356190625.0/4194304.0+(16368226875.0/4194304.0-1588015611.0/4194304.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	case  9 : return    5323807125.0/33554432.0+(-17249135085.0/8388608.0+(212974014825.0/16777216.0+(-398778505125.0/8388608.0+(4004263501875.0/33554432.0+(-886150621125.0/4194304.0+(2269858405815.0/8388608.0+(-1060314480225.0/4194304.0+(5742146284875.0/33554432.0+(-687205715625.0/8388608.0+(441942125625.0/16777216.0+(-42876421497.0/8388608.0+15193976525.0/33554432.0*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2)*r2; break;
	}
    return 0.0;
	}

	//
	// Derivative evaluation (unscaled)
	//
    // Creates 4*Pi* value of the derivative for r2 <= 1

    // Derivative ordering:
    // derivativeList[0] = 0th derivative
    // derivativeList[1] =   x derivative
    // derivativeList[2] =   y derivative
    // derivativeList[3] =   z derivative

    //
    // 0 <=  differentiability order <= 9

    void evaluateDerivatives3D(double r2, double x, double y, double z, std::vector <double>& derivativeList) const
    {
    	int size = derivativeList.size();
    	if(size == 1)
    	{
    	switch(order)
    	{
    	case 2 : derivativeList[0] = evaluation3D_2ndOrder(r2); return;
    	case 4 : derivativeList[0] = evaluation3D_4thOrder(r2); return;
    	case 6 : derivativeList[0] = evaluation3D_6thOrder(r2); return;
    	}}

        double  dVal;
    	if(size >= 4)
    	{
    	switch(order)
    	{
    	case 2 : derivativeList[0] = evaluation3D_2ndOrder(r2); dVal = evaluate2xDq3D_2ndOrder(r2); break;
    	case 4 : derivativeList[0] = evaluation3D_4thOrder(r2); dVal = evaluate2xDq3D_4thOrder(r2); break;
    	case 6 : derivativeList[0] = evaluation3D_6thOrder(r2); dVal = evaluate2xDq3D_6thOrder(r2); break;
    	}
    	derivativeList[1] =  dVal*x;
    	derivativeList[2] =  dVal*y;
    	derivativeList[3] =  dVal*z;
    	}
    }


    double        xPos;    // The position of the mollifier
    double        yPos;
    double        zPos;
    double        laplaceCoeff;
    double      radius;     // The radius of the mollifier
    double    strength;    // The strength of the mollifier

    int           order;    // Order of the mollifier = 2, 4 or 6
    int        exponent;    // The exponent of the mollifer
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
