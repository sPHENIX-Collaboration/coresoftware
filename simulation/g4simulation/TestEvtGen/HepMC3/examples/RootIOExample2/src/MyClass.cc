#include "MyClass.h"

MyClass::MyClass():someint(0), event(0) {}

void MyClass::SetEvent(GenEvent* myevt)
{
    event = myevt;
}

GenEvent* MyClass::GetEvent()
{
    return event;
}

void MyClass::SetInt(int theint)
{
    someint = theint;
}

int MyClass::GetInt()
{
    return someint;
}
