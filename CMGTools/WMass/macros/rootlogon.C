{
  gSystem->Load("libFWCoreFWLite.so");  
  AutoLibraryLoader::enable(); // AutoLibraryLoader::enable() in CMSSW_7_4_X and prior releases
  //FWLiteEnabler::enable();
}
