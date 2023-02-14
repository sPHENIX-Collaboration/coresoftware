/*!
 * \file PHG4MicromegasSurvey.cc
 * \author Hugo Pereira Da Costa <hugo.pereira-da-costa@cea.fr>
 * \brief implements survey data for TPOT definition
 */

#include "PHG4MicromegasSurvey.h"

#include <Geant4/G4ThreeVector.hh>    
#include <Geant4/G4RotationMatrix.hh>    
#include <Geant4/G4SystemOfUnits.hh>

//____________________________________________________________________________________________________
PHG4MicromegasSurvey::PHG4MicromegasSurvey()
{
  // map layer and tile number to detector name
  m_tile_map =  {
    {{55,0},"M5P"},  {{56,0},"M5Z"},  
    {{55,1},"M8P"},  {{56,1},"M8Z"},  
    {{55,2},"M4P"},  {{56,2},"M4Z"},  
    {{55,3},"M10P"}, {{56,3},"M10Z"}, 
    {{55,4},"M9P"},  {{56,4},"M9Z"},  
    {{55,5},"M2P"},  {{56,5},"M2Z"},  
    {{55,6},"M6P"},  {{56,6},"M6Z"},  
    {{55,7},"M7P"},  {{56,7},"M7Z"}   
  };                
  
  // map module name to transformation
  /*
   * these are constructing by mapping the first and last strip positions as given by GEANT4 default implementation of TPOT
   * to the ones measured by the survey crew
   */
  m_transformation_map = {
    { "M10P", 
      { {0.338859, 0.104012, 0.490932}, 
        { 0.999988, 0.00484036, -0.00107067, 
          -0.0048408, 0.999988, -0.000407714, 
          0.00106868, 0.000412892, 1
        } }
    }, 
    { "M10Z", 
      { {0.203866, -0.00685283, 0.125818}, 
        { 0.999989, 0.00473816, -0.000280186, 
          -0.00473828, 0.999989, -0.000413722, 
          0.000278222, 0.000415045, 1
        } }
    }, 
    { "M4P", 
      { {0.775546, 0.281715, 0.603816}, 
        { 0.999924, 0.0120762, 0.00233767, 
          -0.0120718, 0.999925, -0.00186811, 
          -0.00236006, 0.00183975, 1
        } }
    }, 
    { "M4Z", 
      { {0.691254, 0.25956, 0.233509}, 
        { 0.999925, 0.0120751, 0.0021097, 
          -0.0120712, 0.999925, -0.00186256, 
          -0.00213204, 0.00183695, 1
        } }
    }, 
    { "M8P", 
      { {-0.0668271, 0.366721, -0.147311}, 
        { 0.999998, 0.00167415, 0.00100511, 
          -0.00168143, 0.999972, 0.00728541, 
          -0.000992884, -0.00728708, 1
        } }
    }, 
    { "M8Z", 
      { {-0.0845217, 0.333349, -0.512558}, 
        { 0.999996, 0.00167416, 0.00211155, 
          -0.0016895, 0.999972, 0.00728331, 
          -0.0020993, -0.00728685, 1
        } }
    }, 
    { "M5P", 
      { {0.153985, 0.313936, 0.136144}, 
        { 0.999991, 0.00411177, 0.000719481, 
          -0.00411429, 0.999985, 0.00353468, 
          -0.000704936, -0.00353761, 1
        } }
    }, 
    { "M5Z", 
      { {0.0776239, 0.304075, -0.240962}, 
        { 0.999991, 0.00411421, 0.000868929, 
          -0.00411726, 0.999985, 0.00353914, 
          -0.000854355, -0.00354269, 1
        } }
    }, 
    { "M2P", 
      { {-1.16822, 0.745588, 0.702484}, 
        { 0.999958, -0.00915733, 0.000335848, 
          0.00915747, 0.999958, -0.000404121, 
          -0.000332133, 0.00040718, 1
        } }
    }, 
    { "M2Z", 
      { {-1.21001, 0.765289, 0.315568}, 
        { 0.999958, -0.00911263, -0.000571751, 
          0.00911264, 0.999958, 2.19649e-05, 
          0.000571527, -2.71741e-05, 1
        } }
    }, 
    { "M9P", 
      { {-1.2358, 0.82826, 0.350543}, 
        { 0.999936, -0.0110166, 0.00269413, 
          0.0110084, 0.999935, 0.00304423, 
          -0.00272749, -0.00301437, 1
        } }
    }, 
    { "M9Z", 
      { {-1.31984, 0.860159, -0.0225168}, 
        { 0.999936, -0.0110312, 0.00263268, 
          0.0110231, 0.999935, 0.00305265, 
          -0.00266619, -0.00302343, 1
        } }
    }, 
    { "M7P", 
      { {-1.15285, -0.140122, 0.102524}, 
        { 0.999899, -0.0140391, 0.00213771, 
          0.0140361, 0.9999, 0.00143074, 
          -0.00215758, -0.00140059, 1
        } }
    }, 
    { "M7Z", 
      { {-1.20659, -0.241691, -0.269153}, 
        { 0.999899, -0.01402, 0.002198, 
          0.0140168, 0.999901, 0.00146434, 
          -0.00221831, -0.00143339, 1
        } }
    }, 
    { "M6P", 
      { {-0.952177, 0.182855, -0.169552}, 
        { 0.999945, -0.0104378, -0.000665559, 
          0.0104403, 0.999938, 0.00385496, 
          0.00062528, -0.0038617, 1
        } }
    }, 
    { "M6Z", 
      { {-0.992363, -0.0453873, -0.575342}, 
        { 0.999943, -0.0104053, -0.00225987, 
          0.0104119, 0.999942, 0.00292281, 
          0.00222933, -0.00294617, 1
        } }
    }
  };

}

//____________________________________________________________________________________________________
std::string PHG4MicromegasSurvey::get_module_name( int layer, uint tile ) const 
{
  const auto iter = m_tile_map.find( {layer, tile} );
  if( iter == m_tile_map.end() ) 
  {
    std::cout << " PHG4MicromegasSurvey::get_module_name - module not found. layer: " << layer << " tile: " << tile << std::endl;
    return std::string();
  } else return iter->second;
}

//____________________________________________________________________________________________________
G4Transform3D PHG4MicromegasSurvey::get_transformation( int layer, uint tile ) const 
{  
  const auto tile_iter = m_tile_map.find( {layer, tile} );
  if( tile_iter == m_tile_map.end() ) 
  {
    std::cout << " PHG4MicromegasSurvey::get_transformation - module not found. layer: " << layer << " tile: " << tile << std::endl;
    return G4Transform3D();
  }
 
  const auto& tile_name = tile_iter->second;
  const auto transform_iter = m_transformation_map.find(tile_name);
  if( transform_iter == m_transformation_map.end() )
  {
    std::cout 
      << " PHG4MicromegasSurvey::get_transformation - transformation not found."
      << " layer: " << layer 
      << " tile: " << tile
      << " tile_name: " << tile_name 
      << std::endl;
    return G4Transform3D();
  }
  
  // get transformation
  const auto& transformation = transform_iter->second;
  const G4ThreeVector translation( 
    transformation.m_translation[0]*cm, 
    transformation.m_translation[1]*cm, 
    transformation.m_translation[2]*cm);
  const G4RotationMatrix rotation( &transformation.m_rotation[0] );
  return G4Transform3D( rotation, translation );  
}
