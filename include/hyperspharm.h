/**
 * @file hyperspharm.h
 * @author Sylvaus
 * @date Tue July 03 2018
 * @brief
 *
 */

#pragma once

namespace hypershparm
{

class HyperSphericalSurface
{

};

class HyperSphericalCoeffs
{

};

class HyperSpharm
{
  HyperSphericalCoeffs transform(HyperSphericalSurface surface);
  HyperSphericalSurface transform(HyperSphericalCoeffs coeffs);
};

}

