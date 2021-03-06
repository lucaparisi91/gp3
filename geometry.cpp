#include "geometry.h"


size_t geometry::size() const
{
    size_t _size=1;
    for (int i=0;i<AMREX_SPACEDIM;i++)
    {
        _size*=shape[i];
    }

    return _size;
}