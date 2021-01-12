#ifdef __CINT__                            
                                           
#pragma link off all globals;              
#pragma link off all classes;              
#pragma link off all functions;            
#pragma link C++ nestedclasses;            
                                           
#pragma link C++ defined_in "MugastData.h";
#pragma link C++ class std::vector<TVector3>+;
#pragma link C++ class std::vector<std::vector<TVector3>>+;
#pragma link C++ class std::vector<std::vector<TVector3::Theta()>>+;
#pragma link C++ class std::vector<std::vector<double>>+;
#pragma link C++ class MugastData+;
#pragma link C++ class std::vector<MugastData>+;
#endif // __CINT__
