#ifndef GP_FUNCTIONAL_FACTORY_H
#define GP_FUNCTIONAL_FACTORY_H

#include "functional.h"

namespace gp
{
template<class T>
std::shared_ptr<functional> __createFunctional(const json_t & j)
 {
    return std::make_shared<T>(j);
 }

class functionalFactory
{
    public:

    functionalFactory(){}

    
    std::shared_ptr<functional> create(const json_t & j);

    template<class T>
    void registerFunctional(const std::string & key)
    {
        creatorMap[key]= &__createFunctional<T> ; 
    }

    


    private:

    typedef std::shared_ptr<functional> (*functionalCreatorFunc) ( const json_t & j);

    using creatorMap_t = std::map<std::string,functionalCreatorFunc>;
    creatorMap_t creatorMap;

};
}
#endif