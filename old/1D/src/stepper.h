
class eulerStepper
{
public:
  
  template<class op_t>
  void step(state_t state2,state_t state1,Op_t & op,real_t timeStep);
  {
    state2=op(state1)*timeStep;
    state2=stat1  + state2;
  }
  
};
