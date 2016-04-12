// $Id: $                                                                                             
 
/*!
 * \file PROTOTYPE2_FEM.C
 * \brief 
 * \author Jin Huang <jhuang@bnl.gov>
 * \version $Revision:   $
 * \date $Date: $
 */

#include "PROTOTYPE2_FEM.h"

#include <iostream>
#include <cassert>

 int
 PROTOTYPE2_FEM::GetHBDCh(std::string caloname, int i_column, int i_row)
 {
   if (caloname == "HCALIN")
     {
       return 64 + 8 * i_column + 2 * i_row;
     }
   else if (caloname == "HCALOUT")
     {
       return 112 + 8 * i_column + 2 * i_row;
     }
   else if (caloname == "EMCAL")
     {
       // EMcal cable mapping from John haggerty
       assert(i_column >= 0);
       assert(i_column < NCH_EMCAL_COLUMNS);
       assert(i_row >= 0);
       assert(i_row < NCH_EMCAL_ROWS);

       static int canmap[NCH_EMCAL_ROWS*NCH_EMCAL_COLUMNS] =
         {
         // 1 ... 15
             10 + 48, 11 + 48, 8 + 48, 9 + 48, 14 + 48, 15 + 48, 12 + 48, 13
                 + 48,
             // 9 ... 16
             2 + 48, 3 + 48, 0 + 48, 1 + 48, 6 + 48, 7 + 48, 4 + 48, 5 + 48,

             // 17 ... 24
             10 + 32, 11 + 32, 8 + 32, 9 + 32, 14 + 32, 15 + 32, 12 + 32, 13
                 + 32,
             // 25 ... 32
             2 + 32, 3 + 32, 0 + 32, 1 + 32, 6 + 32, 7 + 32, 4 + 32, 5 + 32,

             // 33 ... 40
             10 + 16, 11 + 16, 8 + 16, 9 + 16, 14 + 16, 15 + 16, 12 + 16, 13
                 + 16,
             // 41 42 43 44 45 46 47 48
             2 + 16, 3 + 16, 0 + 16, 1 + 16, 6 + 16, 7 + 16, 4 + 16, 5 + 16,

             // 49 50 51 52 53 54 55 56
             10, 11, 8, 9, 14, 15, 12, 13,
             // 57 58 59 60 61 62 63 64
             2, 3, 0, 1, 6, 7, 4, 5
         };

       const int tower_index = i_column + NCH_EMCAL_COLUMNS * (NCH_EMCAL_ROWS - 1 - i_row);

       assert(tower_index >= 0);
       assert(tower_index < NCH_EMCAL_ROWS*NCH_EMCAL_COLUMNS);

       return canmap[tower_index];
     }

   std::cout << "PROTOTYPE2_FEM::GetHBDCh - invalid input caloname " << caloname
       << " i_column " << i_column << " i_row " << i_row << std::endl;
   exit(1);
   return -9999;
 }
