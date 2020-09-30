#ifndef OP_FACTORY_H
#define OP_FACTORY_H


#include "abstractFactory.h"

/*
Defines a wavefunction factory. Does not anything about jastrows. Id are strings formed by concateneting recursively 'kind' all kind items in the object using '/' as a separator.
*/


class op;

typedef op* (*opCreatorFunc) ( const json_t & j );

template<class func_t>
op * createop(const json_t & j )
{
  return new func_t(j);
}

class opFactory : public abstractFactory<op,std::string, opCreatorFunc>
{
public:
  using abstractFactory_t = abstractFactory<op,std::string, opCreatorFunc>;


  template<class func_t >
  void registerop()
  {
    registerType( func_t::name()  , & (createop<func_t> ) );
  }
 
 op* create(const json_t & j) 
 {
    return abstractFactory_t::create(j["name"].get<std::string>(),j);

 }

};


#endif