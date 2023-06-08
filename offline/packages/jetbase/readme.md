Loopking at JetMapv1.cc::insert,
         specifically with `std::Map<unsigned int, Jet*> type_JetMap`, it makes sense
         to use a TClonesArray here (or even a vector!?!) instead of a map with int -> location

Once inserted, how does the JetMap "own" the Jet memory?
 - I don't think it does. It is just that it delets it under the ::Reset method

 Add parmeters for calculating jet area into FastJetAlgo.h

Store jets in TClonesArray -- have the user supply the array size (if possible!?!), and print a message if it is passed.


Leave in g4jets
