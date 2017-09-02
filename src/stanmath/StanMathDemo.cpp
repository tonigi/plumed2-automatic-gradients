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

#include "StanMathDemo.h"

#define CURVATURE_CITATION "Giorgino T.,  Automatic Gradient Computation Approaches for PLUMED 2  (provisional title)."

using namespace std;

namespace PLMD {
namespace stanmath {

//+PLUMEDOC COLVAR STANMATHDEMO
/*

This file provides a dummy colvar to demonstrate how to implement
functions and have the [Stan Math](https://github.com/stan-dev/math)
library compute their derivatives at compile-time.


\par The Stan Math Library

[Stan Math](https://github.com/stan-dev/math) is a header-only library
providing an extensive library of mathematical and statistical
functions, automatic compile-time differentiation, ODE support, linear
algebra, and so on.

Please see the source code of \ref StanMathDemo.cpp .

> Bob Carpenter, Matthew D. Hoffman, Marcus Brubaker, Daniel Lee,
  Peter Li, and Michael J. Betancourt. 2015. The Stan Math Library:
  Reverse-Mode Automatic Differentiation in C++. arXiv 1509.07164.


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
    // quite readable. [Also: r1 << x(0), x(1), x(2)]
    Eigen::Matrix<T, 3, 1> r1, r2, r3;
    r1 = x.segment(0,3);
    r2 = x.segment(3,3);
    r3 = x.segment(6,3);

    Eigen::Matrix<T, 3, 1> r12, r32, r13;
    r12 = r1-r2;
    r32 = r3-r2;
    r13 = r1-r3;

    T cos2_a = pow(r12.dot(r32),2.0) / r12.dot(r12) / r32.dot(r32);
    T sin2_a = 1.0 - cos2_a;

    T radius2 = r13.dot(r13) / sin2_a / 4.0;
    T radius = sqrt(radius2);

    if(inverse)
      radius = 1.0/radius;

    return(radius);
  }

};


// ----------------------------------------


static const int natoms=3;

class StanMathDemo : public Colvar {
  bool pbc;
  bool inverse;

public:
  static void registerKeywords( Keywords& keys );
  explicit StanMathDemo(const ActionOptions&);
// active methods:
  virtual void calculate();

private:
  void getAtomPositionsAsEigenMatrix(Eigen::Matrix<double,3*natoms,1> &x);
  void setDerivativesFromEigenMatrix(Eigen::Matrix<double,3*natoms,1> x);

};

PLUMED_REGISTER_ACTION(StanMathDemo,"STANMATHDEMO")


// Unpack atom positions as an Eigen matrix object. Passing by
// reference for speed. This function and the next might be made
// static, if we have a way to (a) call the appropriate getPosition()
// method, and (b) know the natoms value for the current colvar. If
// necessary we could use a Dynamic size.
inline void StanMathDemo::getAtomPositionsAsEigenMatrix(Eigen::Matrix<double,3*natoms,1> &x) {
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
inline void StanMathDemo::setDerivativesFromEigenMatrix(Eigen::Matrix<double,3*natoms,1> gx) {
  for (int i=0; i<natoms; i++) {
    int j=i*3;
    Vector v(gx[j],gx[j+1],gx[j+2]);
    setAtomsDerivatives(i,v);
  }
}


void StanMathDemo::registerKeywords( Keywords& keys ) {
  Colvar::registerKeywords( keys );
  keys.add("atoms","ATOMS","the list of three atoms around which to calculate the stanmathdemo");
  keys.addFlag("INVERSE",false,"return the inverse of the radius");
}


StanMathDemo::StanMathDemo(const ActionOptions&ao):
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
void StanMathDemo::calculate() {

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



