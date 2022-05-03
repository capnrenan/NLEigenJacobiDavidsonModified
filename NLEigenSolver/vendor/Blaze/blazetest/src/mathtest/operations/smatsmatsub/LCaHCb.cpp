//=================================================================================================
/*!
//  \file src/mathtest/operations/smatsmatsub/LCaHCb.cpp
//  \brief Source file for the LCaHCb sparse matrix/sparse matrix subtraction math test
//
//  Copyright (C) 2012-2020 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/math/CompressedMatrix.h>
#include <blaze/math/HermitianMatrix.h>
#include <blaze/math/LowerMatrix.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/operations/smatsmatsub/OperationTest.h>
#include <blazetest/system/MathTest.h>

#ifdef BLAZE_USE_HPX_THREADS
#  include <hpx/hpx_main.hpp>
#endif


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'LCaHCb'..." << std::endl;

   using blazetest::mathtest::ScalarA;
   using blazetest::mathtest::ScalarB;

   try
   {
      // Matrix type definitions
      using LCa = blaze::LowerMatrix< blaze::CompressedMatrix<ScalarA> >;
      using HCb = blaze::HermitianMatrix< blaze::CompressedMatrix<ScalarB> >;

      // Creator type definitions
      using CLCa = blazetest::Creator<LCa>;
      using CHCb = blazetest::Creator<HCb>;

      // Running tests with small matrices
      for( size_t i=0UL; i<=6UL; ++i ) {
         for( size_t j=0UL; j<=LCa::maxNonZeros( i ); ++j ) {
            for( size_t k=0UL; k<=i*i; ++k ) {
               RUN_SMATSMATSUB_OPERATION_TEST( CLCa( i, j ), CHCb( i, k ) );
            }
         }
      }

      // Running tests with large matrices
      RUN_SMATSMATSUB_OPERATION_TEST( CLCa(  67UL,  7UL ), CHCb(  67UL, 13UL ) );
      RUN_SMATSMATSUB_OPERATION_TEST( CLCa( 128UL, 16UL ), CHCb( 128UL,  8UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse matrix/sparse matrix subtraction:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************