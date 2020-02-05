
////////////////////////
// USER PARAMETERS
///////////////////////
string __subfolder__ = "figures/reco/";
string __lepton_observable_name__ = "pmu_reco_4mom";
string __nucleon_observable_name__ = "hm_pprot_reco_4mom";
// string __subfolder__ = "figures/true/";
// string __lepton_observable_name__ = "pmu_4mom";
// string __nucleon_observable_name__ = "hm_pprot_4mom";
string __outfile_appendix__ = "";
// string __cuts_string__ = "isRecoChargedLepton == 1 && flagCC0pi == 1 && isRecoProton == 1";
string __cuts_string__ = "isRecoChargedLepton == 1 && isRecoProton == 1 && isRecoPip == 0 && isRecoPim == 0 && Npi0 == 0 &&  && Nother==0";

string __input_data_folder__ = "stephen_output";


/////////////////////////
// GLOBAL VARIABLES
/////////////////////////
vector<string> file_path_list;
vector<string> draw_var_list;
map<string,string> draw_var_to_expr;
map<string, TH1D*> TH1D_buffer;

TLegend* __legend__;
TCanvas* __c_draft__ = new TCanvas("c_draft", "c_draft", 800, 600);
TCanvas* __c_model_comparison__ = new TCanvas("c_model_comparison", "c_model_comparison", 800, 600);
TCanvas* __c_correlation__ = new TCanvas("c_correlation", "c_correlation", 800, 600);

TTree* selectedEvents;
string __input_file_name__;
string __drawing_variable__;

string __draw_string__;

string ERROR    = "\033[1;31m[ERROR] \033[00m";
string INFO  = "\033[1;32m[INFO] \033[00m";
string WARNING   = "\033[1;33m[WARNING] \033[00m";
string ALERT = "\033[1;35m[ALERT] \033[00m";
Color_t colorCells[10] = {kOrange+1, kGreen-3, kTeal+3, kAzure+7, kCyan-2, kBlue-7, kBlue+2, kOrange+9, kRed+2, kPink+9};


/////////////////////////
// FORWARD DECLARATION
/////////////////////////
void simple_plotter(); //main
void do_plot();
void save_canvas(TCanvas *canvas_, string file_name_, string sub_folder_ = __subfolder__);

int __global_plot_line_style__ = 1;













/////////////////////////
// MAIN
/////////////////////////
void simple_plotter(){ // not so simple now...

  file_path_list.emplace_back("neut_lfg_t2k_ch.root");
  file_path_list.emplace_back("neut_rfg_t2k_ch.root");
  file_path_list.emplace_back("neut_sf_t2k_ch.root");
  // file_path_list.emplace_back("neut_RFG_EbLow_t2k_ch.root");
  // file_path_list.emplace_back("neut_RFG_EbHigh_t2k_ch.root");

  draw_var_list.emplace_back("delta_pt");
  draw_var_list.emplace_back("delta_alpha_T");
  draw_var_list.emplace_back("E_nu_QE");
  draw_var_list.emplace_back("E_nu_QE_bias");
  draw_var_list.emplace_back("E_nu_mu_p");
  draw_var_list.emplace_back("E_nu_mu_p_bias");
  draw_var_list.emplace_back("fermi_momentum");

  map<string, TCanvas*> c_model_comparison_list;
  for(auto const draw_var : draw_var_list) {
    string canvas_name = "c_model_comparison_" + draw_var;
    c_model_comparison_list[draw_var] = new TCanvas(canvas_name.c_str(), canvas_name.c_str());
  }

  for(auto const file_path : file_path_list){

    __input_file_name__ = file_path;
    cout << INFO << "> File " << __input_file_name__ << endl;

    int draw_var_count = 0;
    for(auto const draw_var : draw_var_list){

      __drawing_variable__ = draw_var;
      cout << INFO << "  > Variable " << __drawing_variable__ << endl;

      __c_model_comparison__ = c_model_comparison_list[draw_var];
      do_plot();
      draw_var_to_expr[__drawing_variable__] = __draw_string__; // computed in do_plot()

      draw_var_count++;
      if( draw_var_count == draw_var_list.size() ){
        __global_plot_line_style__++;
      }

    }

    TH2D* correlation_plot = new TH2D("correlation_plot", "", 100, 0, 1000, 100, 0, TMath::Pi());
    selectedEvents->Draw(
      (draw_var_to_expr["delta_alpha_T"] + ":" + draw_var_to_expr["delta_pt"] + ">>correlation_plot").c_str(),
      __cuts_string__.c_str(),
      "goff"
    );
    __c_correlation__->cd();
    correlation_plot->SetTitle(__input_file_name__.c_str());
    correlation_plot->GetXaxis()->SetTitle("#deltap_{T}");
    correlation_plot->GetYaxis()->SetTitle("#delta#alpha_{T}");
    correlation_plot->Draw("COLZ");
    string out_name = "correlation_" + __input_file_name__;
    save_canvas(__c_correlation__, out_name);

    // Draw other comparisons
    TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->Scale(1./TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->Integral());
    TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->Scale(1./TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->Integral());
    double peak_value = TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->GetMaximum();
    if(peak_value < TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->GetMaximum()) peak_value = TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->GetMaximum();
    c_model_comparison_list[__input_file_name__ + "_E_nu_bias_comparison"] = new TCanvas(
        (__input_file_name__ + "_E_nu_bias_comparison").c_str(),
        (__input_file_name__ + "_E_nu_bias_comparison").c_str(),
        800, 600);
    c_model_comparison_list[__input_file_name__ + "_E_nu_bias_comparison"]->cd();
    TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->GetYaxis()->SetRangeUser(0, peak_value*1.2);
    TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->Draw("HIST");
    TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->Draw("HIST SAME");

    TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->SetTitle("CCQE Formula");
    TH1D_buffer[__input_file_name__ + "_E_nu_QE_bias_all"]->SetLineColor(colorCells[0]);
    TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->SetTitle("Calorimetric Formula");
    TH1D_buffer[__input_file_name__ + "_E_nu_mu_p_bias_all"]->SetLineColor(colorCells[7]);
    gPad->BuildLegend();
    gPad->SetGridx();
    save_canvas(c_model_comparison_list[__input_file_name__ + "_E_nu_bias_comparison"], __input_file_name__ + "_E_nu_bias_comparison");

  }

  for(auto const draw_var : draw_var_list){
    string canvas_name = "c_model_comparison_" + draw_var;
    c_model_comparison_list[draw_var]->cd();
    gPad->BuildLegend();
    gPad->SetGridx();
    save_canvas(c_model_comparison_list[draw_var], canvas_name);
  }

  __c_draft__->cd();

}












/////////////////////////
// FUNCTIONS
/////////////////////////
void do_plot(){

  gStyle->SetOptStat(0);
  int color_index = 1;

  TFile* file = TFile::Open((__input_data_folder__ + "/" + __input_file_name__).c_str());
  selectedEvents = (TTree*) file->Get("selectedEvents");

  TCanvas* delta_pt = new TCanvas(__drawing_variable__.c_str(), __drawing_variable__.c_str(), 800, 600);

  string delta_pt_x;
  string delta_pt_y;

  double hist_lower_bound;
  double hist_higher_bound;

  string drawing_variable_title;

  if(__drawing_variable__ == "delta_pt"){
    // delta_pt_x = "(pmu_reco_4mom.Px() + hm_pprot_reco_4mom.Px())";
    // delta_pt_y = "(pmu_reco_4mom.Py() + hm_pprot_reco_4mom.Py())";
    delta_pt_x = "(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px())";
    delta_pt_y = "(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py())";
    __draw_string__ = "TMath::Sqrt( TMath::Power(" + delta_pt_x + ",2)" + " + TMath::Power(" + delta_pt_y + ",2) )";

    hist_lower_bound = 0;
    hist_higher_bound = 1000;

    drawing_variable_title = "#deltap_{T} (MeV)";
  }

  if(__drawing_variable__ == "delta_alpha_T"){
    delta_pt_x = "(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px())";
    delta_pt_y = "(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py())";
    string norm_delta_pt = "TMath::Sqrt( TMath::Power(" + delta_pt_x + ",2)" + " + TMath::Power(" + delta_pt_y + ",2) )";
    string norm_pmu = "TMath::Sqrt( TMath::Power(" + __lepton_observable_name__ + ".Px(),2) + TMath::Power(" + __lepton_observable_name__ + ".Py(),2) )";
    string cos_delta_alpha = "-( " + __lepton_observable_name__ + ".Px() * " + delta_pt_x + " + " + __lepton_observable_name__ + ".Py() * " + delta_pt_y + ")/(" + norm_delta_pt + " * " + norm_pmu + ")";
    __draw_string__ = "TMath::ACos(" + cos_delta_alpha + ")";

    hist_lower_bound = 0;
    hist_higher_bound = TMath::Pi();

    drawing_variable_title = "#delta#alpha_{T}";
  }

  if(__drawing_variable__ == "E_nu_QE"){
    string proton_mass = "938.272";
    string neutron_mass = "939.565";
    string muon_mass = "105.658";
    string Eb = "25";

    string numerator = "(";
    numerator += "TMath::Power(" + proton_mass + ",2)";
    numerator += " - ";
    numerator += "TMath::Power(" + muon_mass + ",2)";
    numerator += " + ";
    numerator += "2 * " + __lepton_observable_name__ + ".Energy() * (" + neutron_mass + " - " + Eb + ")";
    numerator += " - ";
    numerator += "TMath::Power(" + neutron_mass + " - " + Eb + ",2)";
    numerator += ")";

    string denominator = "(";
    denominator += "2*(";
    denominator += "(" + neutron_mass + " - " + Eb + ") - " + __lepton_observable_name__ + ".Energy() + " + __lepton_observable_name__ + ".Pz()";
    denominator += ")";
    denominator += ")";

    string E_nu_QE = numerator + "/" + denominator;
    __draw_string__ = E_nu_QE;

    hist_lower_bound = 100;
    hist_higher_bound = 1500;

    drawing_variable_title = "QE Reconstruction of E_{#nu} (MeV)";
  }
  if(__drawing_variable__ == "E_nu_QE_bias"){
    string proton_mass = "938.272";
    string neutron_mass = "939.565";
    string muon_mass = "105.658";
    string Eb = "25";

    string numerator = "(";
    numerator += "TMath::Power(" + proton_mass + ",2)";
    numerator += " - ";
    numerator += "TMath::Power(" + muon_mass + ",2)";
    numerator += " + ";
    numerator += "2 * " + __lepton_observable_name__ + ".Energy() * (" + neutron_mass + " - " + Eb + ")";
    numerator += " - ";
    numerator += "TMath::Power(" + neutron_mass + " - " + Eb + ",2)";
    numerator += ")";

    string denominator = "(";
    denominator += "2*(";
    denominator += "(" + neutron_mass + " - " + Eb + ") - " + __lepton_observable_name__ + ".Energy() + " + __lepton_observable_name__ + ".Pz()";
    denominator += ")";
    denominator += ")";

    string E_nu_QE = numerator + "/" + denominator;
    __draw_string__ = "100*(" + E_nu_QE + " - Enu_true)/Enu_true";

    hist_lower_bound = -50;
    hist_higher_bound = 50;

    drawing_variable_title = "(E_{#nu}^{rec} - E_{#nu}^{true})/E_{#nu}^{true} (%)";
  }
  if(__drawing_variable__ == "E_nu_mu_p"){
    string proton_mass = "938.272";

    string E_nu_mu_p = "(";
    E_nu_mu_p += "" + __lepton_observable_name__ + ".Energy()";
    E_nu_mu_p += " + ";
    E_nu_mu_p += "0.5*( " + __nucleon_observable_name__ + ".Px()*" + __nucleon_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Py()*" + __nucleon_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Pz()*" + __nucleon_observable_name__ + ".Pz() )/" + proton_mass;
    E_nu_mu_p += ")";

    __draw_string__ = E_nu_mu_p;

    hist_lower_bound = 100;
    hist_higher_bound = 1500;

    drawing_variable_title = "Calorimetric Reconstruction of E_{#nu} (MeV)";
  }
  if(__drawing_variable__ == "E_nu_mu_p_bias"){
    string proton_mass = "938.272";

    string E_nu_mu_p = "(";
    E_nu_mu_p += "" + __lepton_observable_name__ + ".Energy()";
    E_nu_mu_p += " + ";
    E_nu_mu_p += "0.5*( " + __nucleon_observable_name__ + ".Px()*" + __nucleon_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Py()*" + __nucleon_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Pz()*" + __nucleon_observable_name__ + ".Pz() )/" + proton_mass;
    E_nu_mu_p += ")";

    __draw_string__ = "100*(" + E_nu_mu_p + " - Enu_true)/Enu_true";

    hist_lower_bound = -50;
    hist_higher_bound = 50;

    drawing_variable_title = "(E_{#nu}^{rec} - E_{#nu}^{true})/E_{#nu}^{true} (%)";
  }
  if(__drawing_variable__ == "fermi_momentum"){

    // string norm_pmu_T = "TMath::Sqrt( TMath::Power(" + __lepton_observable_name__ + ".Px(),2) + TMath::Power(" + __lepton_observable_name__ + ".Py(),2) )";
    string norm_pmu_L = __lepton_observable_name__ + ".Pz()";
    string E_mu = __lepton_observable_name__ + ".Energy()";
    // string norm_pprot_T = "TMath::Sqrt( TMath::Power(" + __nucleon_observable_name__ + ".Px(),2) + TMath::Power(" + __nucleon_observable_name__ + ".Py(),2) )";
    string norm_pprot_L = __nucleon_observable_name__ + ".Pz()";
    string E_prot = __nucleon_observable_name__ + ".Energy()";

    string proton_mass = "938.272";
    string neutron_mass = "939.565";

    // EXAMPLE FOR CARBON
    string binding_energy = "92.16";
    string neutron_separation_energy = "20.3"; // Taking the 1p3/2 subshell -> most populated

    string Mass_A = "(6*" + proton_mass + " + 6*" + neutron_mass + " - " + binding_energy + ")";
    string Mass_A_minus_1 = "( " + Mass_A + " - " + neutron_mass + " + " + neutron_separation_energy + " )";

    // string p_fermi_T = "TMath::Sqrt( TMath::Power(" + norm_pmu_T + ",2) + TMath::Power(" + norm_pprot_T + ",2) )";
    string p_fermi_T = "TMath::Sqrt( ";
    p_fermi_T += "TMath::Power(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px(),2)";
    p_fermi_T += " + TMath::Power(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py(),2)";
    p_fermi_T += " )";
    string X_quantity = "(" + Mass_A + " + " + norm_pmu_L + " + " + norm_pprot_L + " - " + E_mu + " - " + E_prot + ")";

    string p_fermi_L = "( ";
    p_fermi_L += "0.5*" + X_quantity + " - ( TMath::Power(" + p_fermi_T + ",2) + TMath::Power(" + Mass_A_minus_1 + ",2) )/(2*" + X_quantity + ")";
    p_fermi_L += " )";

    string fermi_momentum = "TMath::Sqrt( TMath::Power(" + p_fermi_T + ",2) + TMath::Power( " + p_fermi_L + ",2) )";
    __draw_string__ = fermi_momentum;

    hist_lower_bound = 0;
    hist_higher_bound = 1000;

    drawing_variable_title = "Reconstructed neutron momentum (MeV)";
  }


  int nb_events;

  TH1D* delta_pT_hist = new TH1D("delta_pT_hist", __cuts_string__.c_str(), 100, hist_lower_bound, hist_higher_bound);
  nb_events = selectedEvents->Draw(
    (__draw_string__ + ">>delta_pT_hist").c_str(),
    __cuts_string__.c_str(),
    "goff"
  );
  int norm = delta_pT_hist->GetMaximum();
  delta_pT_hist->Scale(1./norm);
  delta_pT_hist->SetTitle("SuperFGD (all, 100%)");
  delta_pT_hist->SetLineColor(colorCells[color_index]); color_index++;
  delta_pT_hist->SetLineWidth(3);

  double total_nb_events = nb_events;
  delta_pT_hist->GetXaxis()->SetTitle(drawing_variable_title.c_str());
  delta_pT_hist->GetYaxis()->SetTitle("Counts (a.u.)");
  delta_pT_hist->GetYaxis()->SetRangeUser(0, delta_pT_hist->GetMaximum()*1.3);


  TH1D* delta_pT_CCQE_hist = new TH1D("delta_pT_CCQE_hist", __cuts_string__.c_str(), 100, hist_lower_bound, hist_higher_bound);
  nb_events = selectedEvents->Draw(
    (__draw_string__ + ">>delta_pT_CCQE_hist").c_str(),
    (__cuts_string__ + " && Mode == 1").c_str(),
    "goff"
  );
  delta_pT_CCQE_hist->Scale(1./norm);
  delta_pT_CCQE_hist->SetTitle(Form("SuperFGD (CCQE, %d%%)", int((100.*nb_events)/total_nb_events)));
  delta_pT_CCQE_hist->SetLineColor(colorCells[color_index]); color_index++;
  delta_pT_CCQE_hist->SetLineWidth(3);

  TH1D* delta_pT_other_hist = new TH1D("delta_pT_other_hist", __cuts_string__.c_str(), 100, hist_lower_bound, hist_higher_bound);
  nb_events = selectedEvents->Draw(
    (__draw_string__ + ">>delta_pT_other_hist").c_str(),
    (__cuts_string__ + " && Mode > 2").c_str(),
    "goff"
  );
  delta_pT_other_hist->Scale(1./norm);
  delta_pT_other_hist->SetTitle(Form("SuperFGD (Others, %d%%)", int((100.*nb_events)/total_nb_events)));
  delta_pT_other_hist->SetLineColor(colorCells[color_index]); color_index++;
  delta_pT_other_hist->SetLineWidth(3);

  TH1D* delta_pT_2p2h_hist = new TH1D("delta_pT_2p2h_hist", __cuts_string__.c_str(), 100, hist_lower_bound, hist_higher_bound);
  nb_events = selectedEvents->Draw(
    (__draw_string__ + ">>delta_pT_2p2h_hist").c_str(),
    (__cuts_string__ + " && Mode == 2").c_str(),
    "goff"
  );
  delta_pT_2p2h_hist->Scale(1./norm);
  delta_pT_2p2h_hist->SetTitle(Form("SuperFGD (2p2h, %d%%)", int((100.*nb_events)/total_nb_events)));
  delta_pT_2p2h_hist->SetLineColor(colorCells[color_index]); color_index++;
  delta_pT_2p2h_hist->SetLineWidth(3);

  delta_pt->cd();
  delta_pT_hist->Draw("HIST");
  delta_pT_CCQE_hist->Draw("HIST SAME");
  delta_pT_2p2h_hist->Draw("HIST SAME");
  delta_pT_other_hist->Draw("HIST SAME");

  TH1D_buffer[__input_file_name__ + "_" + __drawing_variable__ + "_all"] = (TH1D*) delta_pT_hist->Clone();
  TH1D_buffer[__input_file_name__ + "_" + __drawing_variable__ + "_CCQE"] = (TH1D*) delta_pT_CCQE_hist->Clone();
  TH1D_buffer[__input_file_name__ + "_" + __drawing_variable__ + "_2p2h"] = (TH1D*) delta_pT_2p2h_hist->Clone();
  TH1D_buffer[__input_file_name__ + "_" + __drawing_variable__ + "_other"] = (TH1D*) delta_pT_other_hist->Clone();

  // __legend__ = gPad->BuildLegend(0.55, 0.6, 0.92, 0.92);
  gPad->BuildLegend();
  gPad->SetGridx();
  delta_pT_hist->SetTitle(__input_file_name__.c_str());

  string out_fig_name = __drawing_variable__ + "_" + __input_file_name__;
  save_canvas(delta_pt, out_fig_name);

  __c_model_comparison__->cd();
  delta_pT_hist->SetTitle(__input_file_name__.c_str());
  delta_pT_hist->SetLineStyle(__global_plot_line_style__);
  delta_pT_hist->SetLineColor(colorCells[__global_plot_line_style__]);
  if(__global_plot_line_style__ == 1 ) delta_pT_hist->Draw("HIST");
  else delta_pT_hist->Draw("HIST SAME");

}




void save_canvas(TCanvas *canvas_, string file_name_, string sub_folder_) {

  cout << WARNING << "Saving canvas as " << file_name_ << endl;

  vector<string> extensions;
  extensions.emplace_back("pdf");
  extensions.emplace_back("png");
  extensions.emplace_back("root");
  extensions.emplace_back("C");

  auto old_verbosity = gErrorIgnoreLevel;
  gErrorIgnoreLevel = kFatal;
  for(int i_ext = 0 ; i_ext < int(extensions.size()) ; i_ext++){
    std::system(("mkdir -p ./" + sub_folder_ + "/" + extensions[i_ext]).c_str());
    // cout << WARNING << "Saving as : " << file_name_ << extensions[i_ext] << endl;
    stringstream outpath;
    outpath << "./" << sub_folder_ << "/" << extensions[i_ext] << "/" << file_name_;
    if(not __outfile_appendix__.empty()) outpath << "_" << __outfile_appendix__;
    outpath << "." << extensions[i_ext];
    canvas_->SaveAs(outpath.str().c_str());
  }
  gErrorIgnoreLevel = old_verbosity;

}
