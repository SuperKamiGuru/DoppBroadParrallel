//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4NucleiPropertiesTheoreticalTable.hh 72955 2013-08-14 14:23:14Z gcosmo $
//
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//
// ----------------------------------------------------------------
// Class Description
//   Encapsulates Data from W.D. Myers, W.J. Swiatecki, P. Moller
//   and J.R. Nix, 1. Jan. 1995.
//   Atomic Mass Excess.
// ----------------------------------------------------------------
// Remove "theInstance"  by H.Kurashige (12 Dec. 03)
//

#ifndef G4NucleiPropertiesTheoreticalTable_h
#define G4NucleiPropertiesTheoreticalTable_h 1

//#include "globals.hh"
#include <cmath>
#include <iostream>

class G4NucleiProperties;

class G4NucleiPropertiesTheoreticalTable
{
public:

  // Default constructor
  G4NucleiPropertiesTheoreticalTable(){};

public:

  // Destructor
  ~G4NucleiPropertiesTheoreticalTable() { };

  enum  {nEntries = 8979, shortTableSize = 137};

  // Other Operations
  // all methods are private and can be used only by G4NucleiProperties
  //friend class G4NucleiProperties;


 public: // With Description
  // Operation: GetMassExcess
  static double GetMassExcess(int Z, int A);

  // Operation: GetNuclearMass
  static double GetNuclearMass(int Z, int A);

  // Operation: GetAtomicMass
  static double GetAtomicMass(int Z, int A);

  // Operation: GetBindingEnergy
  static double GetBindingEnergy(int Z, int A);

  // Is the nucleus (Z,A) in table?
  static bool IsInTable(int Z, int A);


public:

  // Operation: GetIndex
  static int GetIndex(int Z, int A);

  static double ElectronicBindingEnergy(int Z);


  // Mass Excess
  static const double AtomicMassExcess[nEntries];

  // Table of Z (number of protons) and A (number of nucleons)
  //        indexArray[0][ ] --> Z
  //        indexArray[1][ ] --> A
  static const int indexArray[2][nEntries];

  // Reduced Table of Z for shorter index search.
  //         The index in this table coincide with Z-1
  //         For each Z value shortTable[Z-1] has the index of
  // the 1st occurrence in the indexArray[][]
  static const int shortTable[shortTableSize];

};


#endif






