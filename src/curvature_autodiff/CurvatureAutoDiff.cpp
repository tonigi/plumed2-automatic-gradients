/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2011-2017 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */

#include "CurvatureAutoDiff.h"

#define CURVATURE_CITATION "T. Giorgino, How to Differentiate Collective Variables in Free Energy Codes: Computer-Algebra Code Generation and Automatic Differentiation, Comp. Phys. Comm.  228 (2018) 258–263, doi:10.1016/j.cpc.2018.02.017"

using namespace std;

namespace PLMD {
namespace curvature_autodiff {

//+PLUMEDOC COLVAR CURVATURE_AUTODIFF
/*

This file provides a dummy colvar to demonstrate how to implement
functions and have the [Stan Math](https://github.com/stan-dev/math)
library compute their derivatives at run-time.


\par The Stan Math Library

[Stan Math](https://github.com/stan-dev/math) is a header-only library
providing an extensive library of mathematical and statistical
functions, automatic compile-time differentiation, ODE support, linear
algebra, and so on.

Please see the source code of \ref CurvatureAutoDiff.cpp .

> Bob Carpenter, Matthew D. Hoffman, Marcus Brubaker, Daniel Lee,
  Peter Li, and Michael J. Betancourt. 2015. The Stan Math Library:
  Reverse-Mode Automatic Differentiation in C++. arXiv 1509.07164.


\par Examples

The following input tells PLUMED to print the curvature at
atoms 1,2,3 and its reciprocal. Usage is similar to \ref CURVATURE_CODEGEN.

\verbatim
c1:  CURVATURE_AUTODIFF ATOMS=1,2,3
c1i: CURVATURE_AUTODIFF ATOMS=1,2,3 INVERSE
PRINT ARG=c1, c1i


*/
//+ENDPLUMEDOC


// ----------------------------------------

/* This struct functor thing may be baffling, but the point is to
   implement operator() assuming an arbitrary number of parameters in
   the argument x . The trick is that instead of scalars you should
   assume to have type T. All of Eigen's linear algebra methods and
   Stan's math functions are available.

   Constant parameters can be "baked-in" the constructor as shown by
   the _inverse_ flag below.

   For readability we could use the following, but let's be explicit
   to avoid confusion with Plumed types:

     using Eigen::Matrix;
     using Eigen::Dynamic;

   When C++14 becomes the standard for PLUMED, we may substitute this
   code with a simpler "lambda" expression (see
   PlumedAutoDiff-Lugano-2019 slides).
*/

struct curvature_fun {
private:
  bool inverse;		// List parameters here
public:
  curvature_fun(bool inverse): inverse(inverse) {}

  template <typename T>
  T operator()(const Eigen::Matrix<T, Eigen::Dynamic, 1>& x)
  const {
    // Reshape for convenience. This will make vector expressions
    // quite readable.

    typedef Eigen::Matrix<T, 3, 1> T3;
    T3 a_, o_, b_;
    a_ << x(0) , x(1) , x(2);
    o_ << x(3) , x(4) , x(5);
    b_ << x(6) , x(7) , x(8);
	
    // "auto xx=..." also works, but it's slower in this case.
    T3 xx = a_-o_;
    T3 yy = b_-o_;
    T3 x_y = a_-b_;

    // Plain Eigen vector operations. Scalars need be of type "T"
    T cos2_a = pow(xx.dot(yy),2.0) / xx.dot(xx) / yy.dot(yy);
    T sin2_a = 1.0 - cos2_a;
    
    T radius2 = x_y.dot(x_y) / sin2_a / 4.0;
    T radius = sqrt(radius2);
	
    if(inverse)
      radius = 1.0/radius;

    return(radius);
  }

};


// ----------------------------------------


static const int natoms=3;

class CurvatureAutoDiff : public Colvar {
  bool pbc;
  bool inverse;

public:
  static void registerKeywords( Keywords& keys );
  explicit CurvatureAutoDiff(const ActionOptions&);
// active methods:
  virtual void calculate();

private:
  void getAtomPositionsAsEigenMatrix(Eigen::Matrix<double,3*natoms,1> &x);
  void setDerivativesFromEigenMatrix(Eigen::Matrix<double,3*natoms,1> x);

};

PLUMED_REGISTER_ACTION(CurvatureAutoDiff,"CURVATURE_AUTODIFF")


// Unpack atom positions as an Eigen matrix object. Passing by
// reference for speed. This function and the next might be made
// static, if we have a way to (a) call the appropriate getPosition()
// method, and (b) know the natoms value for the current colvar. If
// necessary we could use a Dynamic size.
inline void CurvatureAutoDiff::getAtomPositionsAsEigenMatrix(Eigen::Matrix<double,3*natoms,1> &x) {
  int j=0;
  for(int i=0; i<natoms; i++) {
    Vector v=getPosition(i);
    x[j]  =v[0];
    x[j+1]=v[1];
    x[j+2]=v[2];
    j +=3;
  }
}

// Unpack a vector from eigen into a number of Plumed Vectors and set derivatives.
inline void CurvatureAutoDiff::setDerivativesFromEigenMatrix(Eigen::Matrix<double,3*natoms,1> gx) {
  for (int i=0; i<natoms; i++) {
    int j=i*3;
    Vector v(gx[j],gx[j+1],gx[j+2]);
    setAtomsDerivatives(i,v);
  }
}


void CurvatureAutoDiff::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of three atoms around which to calculate the curvature");
  keys.addFlag("INVERSE",false,"return the inverse of the radius");
}


CurvatureAutoDiff::CurvatureAutoDiff(const ActionOptions&ao):
  PLUMED_COLVAR_INIT(ao),
  pbc(true),
  inverse(false)
{
  std::vector<AtomNumber> atoms;
  parseAtomList("ATOMS",atoms);
  if(atoms.size()!=natoms)
    error("Wrong number of atoms");

  parseFlag("INVERSE",inverse);

  bool nopbc=!pbc;
  parseFlag("NOPBC",nopbc);
  pbc=!nopbc;

  checkRead();

  log.printf("  Stan Math test code with INVERSE parameter %s\n",inverse?"TRUE":"FALSE");

  log<<"  Bibliography "
     <<plumed.cite(CURVATURE_CITATION) << "\n";

  addValueWithDerivatives();
  setNotPeriodic();

  requestAtoms(atoms);
}


// calculator
void CurvatureAutoDiff::calculate() {

  if(pbc) makeWhole();

  // Prepare the four arguments.
  // 1. The function to differentiate, declared as above.
  curvature_fun f(inverse);

  // 2. The point at which the function is to be computed, given as a
  // 1D Eigen Matrix.
  Eigen::Matrix<double,3*natoms,1> x;
  getAtomPositionsAsEigenMatrix(x);

  // 3. The variable to receive the return value of the function.
  double f_x;

  // 4. The variable to receive the gradient, as an Eigen Matrix.
  Eigen::Matrix<double,Eigen::Dynamic,1> grad_f_x;

  // 5. The evaluation happens here. The derivative formula have been
  // generated at compile time!
  stan::math::gradient(f, x, f_x, grad_f_x);

  // 6. Store the value and the derivatives.
  setValue(f_x);
  setDerivativesFromEigenMatrix(grad_f_x);

  stan::math::set_zero_all_adjoints();

  setBoxDerivativesNoPbc();

}

}
}



