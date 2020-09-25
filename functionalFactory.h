
#include "abstractFactory.h"

/*
Defines a wavefunction factory. Does not anything about jastrows. Id are strings formed by concateneting recursively 'kind' all kind items in the object using '/' as a separator.
*/


class functional;

typedef functional* (*functionalCreatorFunc) ( const json_t & j );

template<class func_t>
functional * createFunctional(const json_t & j )
{
  return new func_t(j);
}

class functionalFactory : public abstractFactory<functional,std::string, functionalCreatorFunc>
{
public:
  using abstractFactory_t = abstractFactory<functional,std::string, functionalCreatorFunc>;


  template<class func_t >
  void registerFunctional()
  {
    registerType( func_t::name()  , & (createFunctional<func_t> ) );
  }


};