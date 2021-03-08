#include <cassert>
#include "traits.h"
#include "wavefunction.h"

namespace gp
{
    
    namespace normalization
    {
        std::vector<Real> norm(const wavefunction & wave);
        void normalize(wavefunction & wave,const std::vector<Real> & norms);
    }

}