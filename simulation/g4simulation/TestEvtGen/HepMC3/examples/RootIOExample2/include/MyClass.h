#ifndef MYCLASS_H
#define MYCLASS_H

#include "HepMC3/GenEvent.h"
using namespace HepMC3;
/** @class MyClass
 *  @brief Sample class for root I/O test
 */
class MyClass {
public:

    /// @brief Default constructor
    MyClass();

    /// @brief Set HepMC event
    void SetEvent(GenEvent*);

    /// @brief Get HepMC event
    GenEvent* GetEvent();

    /// @brief Set someint
    void SetInt(int);

    /// @brief Get someint
    int GetInt();


private:
    int someint;            ///< Test int
    GenEvent* event; ///< Test event
};

#endif
