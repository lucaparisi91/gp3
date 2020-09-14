#include "gtest/gtest.h"
#include "geometry.h"
#include "kineticEnergy.h"

TEST(geometry,base)
{
  geometrySpherical geo(1,100);
  
  ASSERT_EQ(geo.sizeStep(),1/100.);
  ASSERT_EQ(geo.radius(0), 0.01*0.5);
  ASSERT_NEAR(geo.radius(10), 0.1 + geo.sizeStep() * 0.5 ,1e-7);
  

}

TEST(kineticEnergy,base)
{
  const size_t N = 10000;
  geometrySpherical geo(1.,N);

  state_t state(N);
  for (int i=0;i<N;i++)
    {
      auto r = geo.radius(i);
      state(i)=exp(-r*r);
    }
  
  state_t state2(N);

  auto K=kineticEnergy();
  
  K.evaluate(state2,state,geo);
  
  for (int i=0;i<N;i++)
    {
      EXPECT_NEAR(state2(i) , -2*(1-2*geo.radius(i)*geo.radius(i))*state(i)  ,1e-3);
    }
  
  
}


