void NLReverser(){
	TFile* rtfile = TFile::Open("protodune-nl-MC-fitres.root");
	TGraph* gnlc = (TGraph*)rtfile->Get("NLcorrection");
	TGraph* gresp = (TGraph*)rtfile->Get("Response");
	vector<double> x,y;
	for(int i=1; i<4096; i++){
		double ADCideal = i;
		double nlcorrec = gnlc->Eval(i); // ideal - measurement
		double ADCmeas = ADCideal - nlcorrec;
		x.push_back(ADCideal);
		y.push_back(ADCmeas);
	}


	TGraph* gideal2meas = new TGraph(x.size(), &(x[0]), &(y[0]));
	gideal2meas->GetXaxis()->SetTitle("ideal ADC");
	gideal2meas->GetYaxis()->SetTitle("measured ADC");	
	gideal2meas->SetName("NLideal2meas");
	// gideal2meas->Draw();


	TFile* ofile = new TFile("protodune-nl-MC-NLideal2meas.root","recreate");
	gideal2meas->Write();
	gresp->Write();
	ofile->Close();

	rtfile->Close();
}
