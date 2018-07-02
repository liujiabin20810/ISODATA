//=== File Prolog ================================================================
//      This code was developed by NASA, Goddard Space Flight Center, Code 588
//--- Contents -------------------------------------------------------------------
//      File 'Point.h' contains the definition of Point class in C++.
//
//--- Description ----------------------------------------------------------------
//  This calss was originally implemented for ISOCLUS algorithm but one may use it for
//  other image processing purposes.
//-- Notes:-----------------------------------------------------------------------
//
//-- Development History:--------------------------------------------------------
//   Date             Author                Reference
//   Description
//   
//   March 2002       Nargess Memarsadeghi  NASA GSFC, Code 588
//   Initial implementation
//
//--- DISCLAIMER---------------------------------------------------------------
//
//	This software is provided "as is" without any warranty of any kind, either
//	express, implied, or statutory, including, but not limited to, any
//	warranty that the software will conform to specification, any implied
//	warranties of merchantability, fitness for a particular purpose, and
//	freedom from infringement, and any warranty that the documentation will
//	conform to the program, or any warranty that the software will be error
//	free.
//
//	In no event shall NASA be liable for any damages, including, but not
//	limited to direct, indirect, special or consequential damages, arising out
//	of, resulting from, or in any way connected with this software, whether or
//	not based upon warranty, contract, tort or otherwise, whether or not
//	injury was sustained by persons or property or otherwise, and whether or
//	not loss was sustained from or arose out of the results of, or use of,
//	their software or services provided hereunder.
//--- Warning-------------------------------------------------------------------
//    This software is property of the National Aeronautics and Space
//    Administration.  Unauthorized use or duplication of this software is
//    strictly prohibited.  Authorized users are subject to the following
//    restrictions:
//    *   Neither the author, their corporation, nor NASA is responsible for
//        any consequence of the use of this software.
//    *   The origin of this software must not be misrepresented either by
//        explicit claim or by omission.
//    *   Altered versions of this software must be plainly marked as such.
//    *   This notice may not be removed or altered.
//
//=== End File Prolog=============================================================
#ifndef POINT_H
#define POINT_H

#include <fstream>
class Point
{
public:

	Point();
	Point(int dim, double* p);
	Point* AllocPoint(int dim, double* p);
	Point(const Point&);
	~Point();
	double Norm2DistanceSquared(Point* p) const;
	void print(std::ostream* out) const;
	int getDimension() const;

	double* getPoint()  const;
	void setPoint(double* p);

	Point& operator+(const Point&) const;
	Point& operator+(double) const;
	Point& operator-(double) const;
        Point& operator*(double) const; 
	Point& operator/(double) const; 
	Point& operator=(const Point&) ;
	int operator==(const Point&) const;
	float getCoordinate(int c);
	void  setCoordinate(int c, double value);

private:

	/* a point is a 'dimension'*1 array of type float, with values     */
	/* between 1-255, as far as the image is concerned, the point      */ 	
	/* could have been declared as an 'int*', however since in some	   */
        /* cases we need to calculate means and averages, a floating point */
        /* type is required						   */
	double* point;
	int dimension;

};
#endif

