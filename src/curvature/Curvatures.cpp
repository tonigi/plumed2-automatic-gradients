/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2012-2017 The plumed team
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
#include "multicolvar/MultiColvar.h"
#include "core/PlumedMain.h"
#include "core/ActionRegister.h"
#include "Wrapper.h"


#include <string>
#include <cmath>

using namespace std;
using namespace PLMD::multicolvar;

namespace PLMD {
namespace curvature {

//+PLUMEDOC MCOLVAR CURVATURES
/*
Calculate the local curvature at several triplets of atoms or along a polymer.

You can then calculate functions of the distribution of curvatures
such as the minimum, the number less than a certain quantity and so
on, as per PLUMED2's \ref mcolv functionality. See \ref CURVATURE for
the definition of curvature and details.

The list of places where curvature is to be calculated should
be supplied via multiple ATOMS keywords, each taking a triplet.
Alternatively, the POLYMER keyword can be used to build a list of
curvatures along the triplets of consecutive beads in the list
supplied (see the first two examples below).

The INVERT keyword is used to return the inverse of the curvature
radius.

For the purposes of multicolvar spatial distribution (\ref DENSITY and
related keywords), the center of each multicolvar component is assumed
to be at the geometrical center of the corresponding atom triplet.


\par Examples

The following input tells plumed to calculate the curvature radiuses of the
three circles passing through atoms 1,2 and 3; 2, 3 and 4; and 3, 4 and 5.
The minimum of the three curvatures is used as a colvar and printed.
\verbatim
CURVATURES ATOMS1=1,2,3 ATOMS2=2,3,4 ATOMS3=3,4,5 MIN={BETA=0.1} LABEL=d1
\endverbatim

The following input uses the POLYMER keyword and it is equivalent to
the previous one.

\verbatim
CURVATURES POLYMER=1-5 MIN={BETA=0.1} LABEL=d1
\endverbatim

Count the number of beads in a polymer having a local radius of
curvature less than 3.5:

\verbatim
d1:    CURVATURES POLYMER=1-100 LESS_THAN={RATIONAL R_0=3.5} LOWMEM
PRINT ARG=d1.lessthan
\endverbatim


*/
//+ENDPLUMEDOC


#define POLYMER "POLYMER"

class Curvatures : public MultiColvar {
private:
  bool inverse;
  void readPolymerKeyword( std::vector<AtomNumber>& all_atoms );
public:
  static void registerKeywords( Keywords& keys );
  explicit Curvatures(const ActionOptions&);
// active methods:
  virtual double compute( const unsigned& tindex, AtomValuePack& myatoms ) const ;
/// Returns the number of coordinates of the field
  bool isPeriodic() { return false; }
};

PLUMED_REGISTER_ACTION(Curvatures,"CURVATURES")

void Curvatures::registerKeywords( Keywords& keys ) {
  MultiColvar::registerKeywords( keys );
  keys.add("atoms",POLYMER,"list the beads compounding the polymer");
  keys.addFlag("INVERSE",false,"return the inverse of the radius");
  keys.use("ATOMS"); keys.use("ALT_MIN"); keys.use("LOWEST"); keys.use("HIGHEST");
  keys.use("MEAN"); keys.use("MIN"); keys.use("MAX"); keys.use("LESS_THAN");
  keys.use("MORE_THAN"); keys.use("BETWEEN"); keys.use("HISTOGRAM"); keys.use("MOMENTS");
}

void Curvatures::readPolymerKeyword( std::vector<AtomNumber>& all_atoms ) {
  std::vector<AtomNumber> t;
  const int natoms=3;
  parseAtomList(POLYMER, -1, t);
  if(!t.empty()) {
    // See MultiColvar::readAtomsLikeKeyword
    if( all_atoms.size()>0 )
      error("POLYMER and ATOMS are mutually exclusive.");
    if( t.size() < 3)
      error("POLYMER should have at least 3 atoms");

    ablocks.resize(natoms);
    for(unsigned i=1; i<=t.size()-2; i++) {
      // Are we creating useless atoms?
      ablocks[0].push_back( natoms*(i-1)+0 );
      ablocks[1].push_back( natoms*(i-1)+1 );
      ablocks[2].push_back( natoms*(i-1)+2 );
      all_atoms.push_back( t[i-1] );
      all_atoms.push_back( t[i  ] );
      all_atoms.push_back( t[i+1] );
      atom_lab.push_back( std::pair<unsigned,unsigned>( 0, natoms*(i-1)+0 ) );
      atom_lab.push_back( std::pair<unsigned,unsigned>( 0, natoms*(i-1)+1 ) );
      atom_lab.push_back( std::pair<unsigned,unsigned>( 0, natoms*(i-1)+2 ) );
      log.printf("  Polymer bead %d is calculated from curvature of atoms: %d %d %d\n",
                 i-1, t[i-1].serial(), t[i].serial(), t[i+1].serial());
    }
    if( all_atoms.size()>0 ) {
      nblock=0;
      for(unsigned i=0; i<ablocks[0].size(); ++i) addTaskToList( i );
    }
  }
}


Curvatures::Curvatures(const ActionOptions&ao):
  PLUMED_MULTICOLVAR_INIT(ao),
  inverse(false)
{
  // Read in the atoms
  std::vector<AtomNumber> all_atoms;

  int natoms=3;

  readAtomsLikeKeyword( "ATOMS", natoms, all_atoms );
  // was: readAtoms( natoms, all_atoms ); which is not ok because it calls setupMultiColvarBase

  readPolymerKeyword(all_atoms);

  setupMultiColvarBase( all_atoms );

  // Invert flag
  parseFlag("INVERSE",inverse);

  // Using the center of the three atoms as the central atoms. This is
  // consistent with the fact that the colvar is invariant by
  // permutations of the input atoms, but may be somewhat at variance
  // from those assuming the list is ordered along the polymer.
  std::vector<bool> catom_ind(natoms, true);
  setAtomsForCentralAtom( catom_ind );

  // Read in the vessels
  readVesselKeywords();

  // And check everything has been read in correctly
  checkRead();

  log<<"  Bibliography "
     <<plumed.cite(CURVATURE_CITATION)
     <<"\n";

}


double Curvatures::compute( const unsigned& tindex, AtomValuePack& myatoms ) const {

  // PBC?

  Vector a=myatoms.getPosition(0);
  Vector b=myatoms.getPosition(1);
  Vector c=myatoms.getPosition(2);

  double r=curvature_radius(a[0],a[1],a[2],
                            b[0],b[1],b[2],
                            c[0],c[1],c[2]);

  double ga_d[3], gb_d[3], gc_d[3];
  curvature_radius_grad(a[0],a[1],a[2],
                        b[0],b[1],b[2],
                        c[0],c[1],c[2],
                        ga_d, gb_d, gc_d );

  Vector ga(ga_d[0], ga_d[1], ga_d[2]);
  Vector gb(gb_d[0], gb_d[1], gb_d[2]);
  Vector gc(gc_d[0], gc_d[1], gc_d[2]);

  double value;

  if(!inverse) {
    // We can't do much if r=inf
    value=r;
  } else {
    value=1/r;
    if (value > 0.0) {
      double minus_inv_r2 = -1.0/(r*r);
      ga = minus_inv_r2*ga;
      gb = minus_inv_r2*gb;
      gc = minus_inv_r2*gc;
    } else {
      log.printf("CURVATURE: radius %f occurred, setting null gradient\n",r);
      Vector v0(0,0,0);
      ga = gb = gc = v0;
    }
  }

  addAtomDerivatives(1, 0, ga, myatoms);
  addAtomDerivatives(1, 1, gb, myatoms);
  addAtomDerivatives(1, 2, gc, myatoms);

  myatoms.addBoxDerivatives( 1, -(Tensor(a,ga)+Tensor(b,gb)+Tensor(c,gc)) );

  return value;
}

}
}

