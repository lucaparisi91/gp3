#include "functionalFactory.h"

namespace gp
{


std::shared_ptr<functional> functionalFactory::create(const json_t & j)
{
     creatorMap_t::const_iterator i;
     std::string key=j["kind"];
    i=creatorMap.find(key);


        if (i!= creatorMap.end())
        {
	    	auto func=(i->second)( j  );
            
            return func;
        }
         else
        {
	    throw factoryIdNotRecorded(key);
        }

}

}