#include "gpDriver.h"
#include <iostream>
#include "initializer.h"

int main(int argc,char**argv)
{
    if (argc !=2)
    {
        std::cerr << "Input filename not provided" << std::endl;
        exit(1);
    }

    std::string filename(argv[1]);
    

    gp::json_t j;

    std::ifstream f;

    f.open(filename);

    f >> j;
    f.close();

    initializer::getInstance().init();
    gp::gpDriver::run(j);

}