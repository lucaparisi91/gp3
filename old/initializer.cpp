#include "initializer.h"
#include "functional.h"
#include "operators.h"
#include "gpDropletFunctional3D.h"
initializer *initializer::s_instance = 0;

initializer::initializer() : isInitialized(false)
{
    func.registerFunctional<harmonicFunctional>();
    func.registerFunctional<gpDroplet3D>();
    
    // operators
    
    operatorFac.registerop<stencilLaplacianOperator>();
    operatorFac.registerop<amrexLaplacianOperator>();
    operatorFac.registerop<stencilLaplacian<2> >();

}