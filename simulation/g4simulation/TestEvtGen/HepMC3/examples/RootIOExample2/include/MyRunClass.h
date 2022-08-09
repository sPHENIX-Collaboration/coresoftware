#ifndef MYRUNCLASS_H
#define MYRUNCLASS_H

#include "HepMC3/GenRunInfo.h"
using namespace HepMC3;
/** @class MyRunClass
 *  @brief Sample class for root I/O test
 */
class MyRunClass {
public:

    /// @brief Default constructor
    MyRunClass();

    /// @brief Set HepMC event
    void SetRunInfo(GenRunInfo*);

    /// @brief Get HepMC event
    GenRunInfo* GetRunInfo();

    /// @brief Set someint
    void SetInt(int);

    /// @brief Get someint
    int GetInt();

private:
    int someint;            ///< Test int
    GenRunInfo* run; ///< Test run info
};

#endif
