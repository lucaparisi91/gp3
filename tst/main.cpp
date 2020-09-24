#include "gtest/gtest.h"
#include <AMReX_MultiFabUtil.H>

int main(int argc, char **argv) {

  amrex::Initialize(MPI_COMM_WORLD);

  ::testing::TestEventListeners& listeners =
    ::testing::UnitTest::GetInstance()->listeners();
  //if (pTools::rank() != 0) {
  //  delete listeners.Release(listeners.default_result_printer());
  //}
    ::testing::InitGoogleTest(&argc, argv);
    
  int ret= RUN_ALL_TESTS();
  amrex::Finalize();
}


