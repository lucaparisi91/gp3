#include "initializer.h"
#include "functional.h"

initializer *initializer::s_instance = 0;

initializer::initializer() : isInitialized(false)
{
    func.registerFunctional<harmonicFunctional>();

}