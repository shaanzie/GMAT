//$Id$
//------------------------------------------------------------------------------
//                                  SpiceAxes
//------------------------------------------------------------------------------
// GMAT: General Mission Analysis Tool.
//
// Copyright (c) 2002-2025 United States Government as represented by the
// Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License"); 
// You may not use this file except in compliance with the License. 
// You may obtain a copy of the License at:
// http://www.apache.org/licenses/LICENSE-2.0 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either 
// express or implied.   See the License for the specific language
// governing permissions and limitations under the License.
//
// Developed jointly by NASA/GSFC and Pearl River Technologies under
// FDSS III.
//
// Author: Peter Candell
// Created: 2024/04/18
//
/**
 * Implementation of the SpiceAxes class.
 *
 */
 //------------------------------------------------------------------------------

#include <iostream>
#include "gmatdefs.hpp"
#include "GmatBase.hpp"
#include "SpiceAxes.hpp"
#include "AxisSystem.hpp"
#include "SolarSystem.hpp"
#include "CelestialBody.hpp"
#include "Planet.hpp"
#include "RealUtilities.hpp"
#include "Linear.hpp"
#include "GmatConstants.hpp"
#include "Rvector3.hpp"
#include "TimeSystemConverter.hpp"
#include "CoordinateSystemException.hpp"
#include "MessageInterface.hpp"
#include "Attitude.hpp"
#include "Spacecraft.hpp"
#include "AttitudeConversionUtility.hpp"

using namespace GmatMathUtil;        // for trig functions, etc.
using namespace GmatTimeConstants;   // for JD offsets, etc.

//------------------------------------------------------------------------------
// static data
//------------------------------------------------------------------------------
// none

//------------------------------------------------------------------------------
// public methods
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  SpiceAxes(const std::string &itsName);
//------------------------------------------------------------------------------
/**
 * Constructs base SpiceAxes structures
 * (default constructor).
 *
 * @param itsName Optional name for the object.  Defaults to "".
 */
 //------------------------------------------------------------------------------
SpiceAxes::SpiceAxes(const std::string &itsName) :
   AxisSystem("Spice", itsName),
   prevEpoch(0.0),
   prevUpdateInterval(-99.9),
   prevOriginUpdateInterval(-99.9),
   spiceFrameId("")
{
   objectTypeNames.push_back("SpiceAxes");
   parameterCount = SpiceAxesParamCount;

   needsCBOrigin = false; // true;
}


//------------------------------------------------------------------------------
//  SpiceAxes(const SpiceAxes &sAxes);
//------------------------------------------------------------------------------
/**
 * Constructs base SpiceAxes structures used in derived classes, by copying
 * the input instance (copy constructor).
 *
 * @param sAxes  SpiceAxes instance to copy to create "this" instance.
 */
 //------------------------------------------------------------------------------
SpiceAxes::SpiceAxes(const SpiceAxes &sAxes) :
   AxisSystem(sAxes),
   prevEpoch(sAxes.prevEpoch),
   prevUpdateInterval(sAxes.prevUpdateInterval),
   prevOriginUpdateInterval(sAxes.prevOriginUpdateInterval),
   spiceFrameId(sAxes.spiceFrameId)
{
#ifdef DEBUG_BF_RECOMPUTE
   MessageInterface::ShowMessage("Constructing a new SpiceAxes (%p) from the old one (%p)\n",
      this, &sAxes);
   MessageInterface::ShowMessage("   and prevEpoch(old) %12.10f copied to prevEpoch(new) %12.10f\n",
      sAxes.prevEpoch.GetMjd(), prevEpoch.GetMjd());
#endif
}

//------------------------------------------------------------------------------
//  SpiceAxes& operator=(const SpiceAxes &sAxes)
//------------------------------------------------------------------------------
/**
 * Assignment operator for SpiceAxes structures.
 *
 * @param sAxes The original that is being copied.
 *
 * @return Reference to this object
 */
 //------------------------------------------------------------------------------
const SpiceAxes& SpiceAxes::operator=(const SpiceAxes &sAxes)
{
   if (&sAxes == this)
      return *this;
   AxisSystem::operator=(sAxes);
   prevEpoch = sAxes.prevEpoch;
   prevUpdateInterval = sAxes.prevUpdateInterval;
   prevOriginUpdateInterval = sAxes.prevOriginUpdateInterval;
   spiceFrameId = sAxes.spiceFrameId;

   return *this;
}


//------------------------------------------------------------------------------
//  ~SpiceAxes()
//------------------------------------------------------------------------------
/**
 * Destructor.
 */
 //------------------------------------------------------------------------------
SpiceAxes::~SpiceAxes()
{
}

//------------------------------------------------------------------------------
//  GmatCoordinate::ParameterUsage UsesEopFile(
//                                 const std::string &forBaseSystem) const
//------------------------------------------------------------------------------
/**
 * @see AxisSystem
 */
 //------------------------------------------------------------------------------
GmatCoordinate::ParameterUsage SpiceAxes::UsesEopFile(
   const std::string &forBaseSystem) const
{
   return GmatCoordinate::NOT_USED;
}

//------------------------------------------------------------------------------
//  GmatCoordinate::ParameterUsage UsesItrfFile() const
//------------------------------------------------------------------------------
/**
 * @see AxisSystem
 */
 //------------------------------------------------------------------------------
GmatCoordinate::ParameterUsage SpiceAxes::UsesItrfFile() const
{
   return GmatCoordinate::NOT_USED;
}

//------------------------------------------------------------------------------
//  GmatCoordinate::ParameterUsage UsesNutationUpdateInterval() const
//------------------------------------------------------------------------------
/**
 * @see AxisSystem
 */
 //---------------------------------------------------------------------------
GmatCoordinate::ParameterUsage SpiceAxes::UsesNutationUpdateInterval() const
{
   return GmatCoordinate::NOT_USED;
}


//------------------------------------------------------------------------------
//  bool Initialize()
//------------------------------------------------------------------------------
/**
 * Initialization method for this SpiceAxes.
 *
 * @return success flag
 */
 //------------------------------------------------------------------------------
bool SpiceAxes::Initialize()
{
#ifdef DEBUG_SpiceAxes_INIT
   MessageInterface::ShowMessage("Entering SpiceAxes::Init ...\n");
   if (origin)
      MessageInterface::ShowMessage("origin is %s<%p>\n",
         origin->GetName().c_str(), origin);
#endif
   AxisSystem::Initialize();
   if (originName == GmatSolarSystemDefaults::EARTH_NAME) InitializeFK5();

   // Check for spacecraft origin whose attitude does not compute rates
   if (origin && origin->IsOfType("Spacecraft"))
   {
      Attitude *att = (Attitude*)origin->GetRefObject(Gmat::ATTITUDE, "");
      if (!att)
      {
         std::string errmsg = "The value of \"";
         errmsg += origin->GetName() + "\" for field \"Origin\"";
         errmsg += " on Spice coordinate system \"" + coordName;
         errmsg += "\" is not an allowed value.\n";
         errmsg += "The allowed values are: ";
         errmsg += "[ Spacecraft with an attitude model ].\n";
         throw CoordinateSystemException(errmsg);
      }
      if ((!allowNoRates) && !(att->ModelComputesRates()))
      {
         std::string errmsg = "The value of \"";
         errmsg += origin->GetName() + "\" for field \"Origin\"";
         errmsg += " on Spice coordinate system \"" + coordName;
         errmsg += "\" is not an allowed value.\n";
         errmsg += "The allowed values are: ";
         errmsg += "[ Spacecraft whose attitude model computes rates ].\n";
         throw CoordinateSystemException(errmsg);
      }
   }

#ifdef DEBUG_FIRST_CALL
   firstCallFired = false;
#endif

   return true;
}


//------------------------------------------------------------------------------
// public methods inherited from GmatBase
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//  GmatBase* Clone() const
//------------------------------------------------------------------------------
/**
 * This method returns a clone of the SpiceAxes.
 *
 * @return clone of the SpiceAxes.
 *
 */
 //------------------------------------------------------------------------------
GmatBase* SpiceAxes::Clone() const
{
   return (new SpiceAxes(*this));
}

//------------------------------------------------------------------------------
//  bool SetRefObject(GmatBase *obj, const UnsignedInt type,
//                    const std::string &name)
//------------------------------------------------------------------------------
/**
 * This method sets a reference object for the SpiceAxes class.  This is
 * overridden from the CoordinateBase version, in order to make sure the origin
 * is a CelestialBody.
 *
 * @param obj   pointer to the reference object
 * @param type  type of the reference object
 * @param name  name of the reference object
 *
 * @return true if successful; otherwise, false.
 *
 */
 //------------------------------------------------------------------------------
bool SpiceAxes::SetRefObject(GmatBase *obj, const UnsignedInt type,
   const std::string &name)
{
#ifdef DEBUG_SET_REF
   MessageInterface::ShowMessage
   ("SpiceAxes::SetRefObject() <%s>, obj=<%p>'%s', name=%s\n", GetName().c_str(),
      obj, obj ? obj->GetName().c_str() : "NULL", name.c_str());
#endif

   if (obj == NULL)
      return false;

   if (name == originName)
   {
      if ((!obj->IsOfType("CelestialBody")) && (!obj->IsOfType("Spacecraft")))
      {
         CoordinateSystemException cse("");
         cse.SetDetails(errorMessageFormatUnnamed.c_str(),
            (obj->GetName()).c_str(),
            "Origin", "Celestial Body or Spacecraft");
         throw cse;
      }
   }
#ifdef DEBUG_SET_REF
   MessageInterface::ShowMessage
   ("SpiceAxes::SetRefObject() returning GmatBase::SetRefObject()\n");
#endif
   return AxisSystem::SetRefObject(obj, type, name);
}


//------------------------------------------------------------------------------
//  void SpiceAxes::SetSpiceFrameId(std::string spiceId)
//------------------------------------------------------------------------------
/**
 * Sets the spice frame id which turns on using Spice for rotations.
 * To turn off: set to ""
 *
 * @param spiceId   name of the spice frame to use
 *
 */
 //------------------------------------------------------------------------------
void SpiceAxes::SetSpiceFrameId(std::string spiceId)
{
   spiceFrameId = spiceId;
}

//------------------------------------------------------------------------------
// protected methods
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//  void CalculateRotationMatrix(const A1Mjd &atEpoch, 
//                               bool forceComputation = false)
//------------------------------------------------------------------------------
/**
 * This method will compute the rotMatrix and rotDotMatrix used for rotations
 * from/to this AxisSystem to/from the MJ2000Eq system
 *
 * @param atEpoch          epoch at which to compute the rotation matrix
 * @param forceComputation force computation even if it is not time to do it
 *                         (default is false)
 */
 //------------------------------------------------------------------------------
void SpiceAxes::CalculateRotationMatrix(const A1Mjd &atEpoch,
   bool forceComputation)
{
   CalculateRotationMatrix(GmatTime(atEpoch.Get()), forceComputation);
}

//------------------------------------------------------------------------------
//  void CalculateRotationMatrix(const GmatTime &atEpoch, 
//                               bool forceComputation = false)
//------------------------------------------------------------------------------
/**
 * This method will compute the rotMatrix and rotDotMatrix used for rotations
 * from/to this AxisSystem to/from the MJ2000Eq system
 *
 * @param atEpoch          epoch at which to compute the rotation matrix
 * @param forceComputation force computation even if it is not time to do it
 *                         (default is false)
 */
 //------------------------------------------------------------------------------
void SpiceAxes::CalculateRotationMatrix(const GmatTime &atEpoch,
   bool forceComputation)
{
#ifdef DEBUG_FIRST_CALL
   if (!firstCallFired)
      MessageInterface::ShowMessage(
         "Calling SpiceAxes::CalculateRotationMatrix at epoch %.15lf; \n", atEpoch.GetMjd());
#endif
   Real theEpoch = atEpoch.GetMjd();
#ifdef DEBUG_BF_EPOCHS
   MessageInterface::ShowMessage("SpiceAxes::CalculateRotationMatrix(%s)   epoch = %12.10f, prevEpoch = %12.10f ...... ",
      (coordName.c_str()), theEpoch, prevEpoch.GetMjd());
#endif
#ifdef DEBUG_BF_RECOMPUTE
   MessageInterface::ShowMessage("Entering SpiceAxes::CalculateRotationMatrix on object %s (%p) of type %s, origin = %s\n",
      (coordName.c_str()), this, (GetTypeName()).c_str(), originName.c_str());
   MessageInterface::ShowMessage("     epoch = %12.10f, prevEpoch = %12.10f\n", theEpoch, prevEpoch.GetMjd());
#endif

   // if the origin is a Spacecraft, we need to get attitude data from the
   // spacecraft in order to compute the R and Rdot matrices
   if (origin->IsOfType("Spacecraft"))
   {
      if ((!forceComputation) && (IsEqual(atEpoch, prevEpoch)))
      {
#ifdef DEBUG_BF_RECOMPUTE
         MessageInterface::ShowMessage(
            "Don't need to recompute R and Rdot SpiceAxes for Spacecraft "
            "%s at this time!!\n",
            origin->GetName().c_str());
#endif
         return;
      }
#ifdef DEBUG_BF_SC
      if (!firstCallFired)
      {
         MessageInterface::ShowMessage(
            "ABOUT to get dcm (and possibly av) from spacecraft %s\n",
            origin->GetName().c_str());
      }
#endif
      // First compute the R matrix
      Spacecraft *sc = (Spacecraft*)origin;
      // DCM from spacecraft is inertial-to-body
      Rmatrix33  dcm = sc->GetAttitude(theEpoch);

      // we want the rotation matrix to be body-to-inertial
      rotMatrix = dcm.Transpose();

      Attitude *att = (Attitude*)sc->GetRefObject(Gmat::ATTITUDE, "");
      if (!att->ModelComputesRates())
      {
         if (allowNoRates)
         {
            rotDotMatrix = Rmatrix33();  // zero matrix
         }
         else
         {
            std::string errmsg = "The value of \"";
            errmsg += origin->GetName() + "\" for field \"Origin\"";
            errmsg += " on Spice coordinate system \"" + coordName;
            errmsg += "\" is not an allowed value.\n";
            errmsg += "The allowed values are: ";
            errmsg += "[ Spacecraft whose attitude model computes rates ].\n";
            throw CoordinateSystemException(errmsg);
         }
      }
      else
      {
         // av from spacecraft is wrt inertial, expressed in inertial
         Rvector3   av = sc->GetAngularVelocity(theEpoch);
#ifdef DEBUG_BF_SC
         if (!firstCallFired)
         {
            MessageInterface::ShowMessage("Computing R and Rdot for SpiceAxes with spacecraft origin\n");
            MessageInterface::ShowMessage("   Spacecraft = %s\n", sc->GetName().c_str());
            MessageInterface::ShowMessage("   dcm    = %s\n", dcm.ToString().c_str());
            MessageInterface::ShowMessage("   av     = %s\n", av.ToString().c_str());
         }
#endif

         // rotate angular velocity from inertial to body
         Rvector3  avB = dcm * av;
         // compute Rdot
         Rmatrix33 skew(0.0, -avB[2], avB[1],
            avB[2], 0.0, -avB[0],
            -avB[1], avB[0], 0.0);
         Rmatrix33 RdotBI = -skew * dcm;
         // we want Rdot to be body-to-inertial
         rotDotMatrix = RdotBI.Transpose();
#ifdef DEBUG_BF_SC
         if (!firstCallFired)
         {
            MessageInterface::ShowMessage("   avB    = %s\n", avB.ToString().c_str());
            MessageInterface::ShowMessage("   skew   = %s\n", skew.ToString().c_str());
            MessageInterface::ShowMessage("   RdotBI = %s\n", RdotBI.ToString().c_str());
         }
#endif
      }
#ifdef DEBUG_BF_SC
      if (!firstCallFired)
      {
         MessageInterface::ShowMessage("   rotMatrix    = %s\n", RdotBI.ToString().c_str());
         MessageInterface::ShowMessage("   rotDotMatrix = %s\n", RdotBI.ToString().c_str());
      }
#endif
   }
   else
   {
      // Code must check to see if we need to recompute.  Recomputation is only necessary
      // if one or more of the following conditions are true:
      // 1. the epoch is different (within machine precision) from the epoch at the last computation
      // 2. if the origin is the Earth, the nutation update interval has changed
      // 3. if the origin is Luna, the rotation data source has changed

      if ((!forceComputation) &&
         (originName == GmatSolarSystemDefaults::MOON_NAME) &&
         (IsEqual(atEpoch, prevEpoch)))
      {
#ifdef DEBUG_BF_RECOMPUTE
         MessageInterface::ShowMessage("Don't need to recompute for Luna at this time!!\n");
#endif
         return;
      }

      // compute rotMatrix and rotDotMatrix
      if (spiceFrameId != "" || ((CelestialBody*)origin)->GetRotationDataSource() == Gmat::SPICE_KERNEL)
      {
         std::string spiceToFrame;

         if (spiceFrameId != "")
            spiceToFrame = spiceFrameId;
         else //((CelestialBody*)origin)->GetRotationDataSource() == Gmat::SPICE_KERNEL)
            spiceToFrame = origin->GetStringParameter("SpiceFrameId");
         CalculateSpiceFrameRotationMatrix(spiceToFrame, atEpoch);
         
      }
      else
      {
         std::string errmsg = "Attempting use a Spice frame to calculate a rotation matrix but no Spice frame was specified for either oirgin: \"";
         errmsg += origin->GetName();
         errmsg += "\" or coordinate system \"" + coordName + "\"";
         throw CoordinateSystemException(errmsg);
      }
   }
}