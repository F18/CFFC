  string mech_name = "CH4_WD1b";
  string mech_file = "CH4.xml";
  double Sc[5] = {1.0, 1.0, 1.0, 1.0, 1.0};
  double y[5] = {0.2, 0.2, 0.2, 0.2, 0.2};
  bool constant_schmidt = true;
  Mixture::setMixture(mech_name,
		      mech_file,
		      Sc, 
		      constant_schmidt);
  Mixture M;
  M.setState_TPY(300.0, 101325.0, y);

  cout << "HERE\n";
  cout << M;
  cout << "HERE\n";

  Flame2D_pState::setMixture(mech_name,
			     mech_file,
			     Sc, 
			     constant_schmidt);
  Flame2D_pState W;

  cout << "HERE\n";
  cout << W;
  cout << "HERE\n";

