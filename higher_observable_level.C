#include<dirent.h>

string __input_data_folder__ = "stephen_output";
string __output_folder__ = "output";


TFile* __current_tfile__;
TFile* __outfile_tfile__;
TTree* __current_ttree__;
TTree* __outfile_ttree__;
TEntryList* __entry_list__;
string __cuts_string__;
map<string, string> __variable_name_to_formula_map__;
vector<string> __variables_list__;
vector<string> __expressions_list__;
map<string, double> __observable_values__;
vector<TTreeFormula*> __formula_list__;

string ERROR    = "\033[1;31m[ERROR] \033[00m";
string INFO  = "\033[1;32m[INFO] \033[00m";
string WARNING   = "\033[1;33m[WARNING] \033[00m";
string ALERT = "\033[1;35m[ALERT] \033[00m";

// main
void higher_observable_level();
void fill_variables_expression();
void hook_tformula_on_input_tree();

// tools
int verbosity_level = 0;
double __last_displayed_value__ = -1;
vector<std::string> get_list_of_files_in_directory(std::string folder_path_="./", std::string file_format_="*");
bool do_string_has_valid_format(std::string string_to_test_, std::string string_format_);
std::string join_vector_string(std::vector<std::string> string_list_, std::string delimiter_, int begin_index_ = 0 , int end_index_ = 0);
std::vector<std::string> split_string(std::string input_string_, char delimiter_);
std::vector<std::string> split_string(std::string input_string_, std::string delimiter_);
bool do_tfile_is_valid(TFile *input_tfile_, bool check_if_writable_ = false);
void display_loading(int current_index_, int end_index_, string title_ = "Processing", bool force_display_ = false);


void higher_observable_level(){

  cout << INFO << "Input folder : " << __input_data_folder__ << endl;
  cout << INFO << "Output folder : " << __output_folder__ << endl;

  // vector<string> files_list = get_list_of_files_in_directory(__input_data_folder__, "*.root");
  vector<string> files_list;
  files_list.emplace_back("neut_rfg_t2k_ch.root");
  fill_variables_expression();

  for( auto const &file : files_list ){

    cout << INFO << "Processing : " << file << endl;

    __current_tfile__ = TFile::Open((__input_data_folder__ + "/" + file).c_str());
    do_tfile_is_valid(__current_tfile__);
    __current_ttree__ = (TTree*) __current_tfile__->Get("selectedEvents");

    hook_tformula_on_input_tree();

    std::system(("mkdir -p ./" + __output_folder__).c_str());
    __outfile_tfile__ = TFile::Open((__output_folder__ + "/HO_" + file).c_str(), "RECREATE");
    __outfile_tfile__->cd();
    __outfile_ttree__ = new TTree("selectedEvents", "selectedEvents");

    // Initializing the data handlers (attributing a memory location for all variables)
    for(auto const &variable_name : __variables_list__){
      __observable_values__[variable_name] = 0.;
    }

    // Hook values to the output ttrees
    vector<string> sorted_variables_name_list(__variables_list__);
    std::sort(sorted_variables_name_list.begin(), sorted_variables_name_list.end());
    string variable_name = "";
    cout << WARNING << "The following variables will be processed :" << endl;
    for(int i_var = 0 ; i_var < int(sorted_variables_name_list.size()) ; i_var++){
      variable_name = sorted_variables_name_list[i_var];
      __outfile_ttree__->Branch(variable_name.c_str(), &__observable_values__[variable_name]);
      cout << WARNING << " " << variable_name << endl;
    }

    // __outfile_ttree__ = __current_ttree__->CloneTree(0);
    TTreeFormula *cuts_formula = new TTreeFormula("__cuts_string__", __cuts_string__.c_str(), __current_ttree__);
    __current_ttree__->SetNotify(cuts_formula); // This is needed only for TChain.
    for (Int_t i=0;i<__current_ttree__->GetEntries(); i++) {
      display_loading(i,__current_ttree__->GetEntries(), WARNING + "Looping over the TTree...");
      __current_ttree__->GetEntry(i); // Does not read any data yet (TTreeFormula will read what it needs).
      Bool_t do_entry_passes_cuts = false;
      for(Int_t j = 0; j < cuts_formula->GetNdata(); ++j) {
        if ( cuts_formula->EvalInstance(j) ) {
          do_entry_passes_cuts = true;
          break;
        }
      }
      if (do_entry_passes_cuts) {
        // eval all variables
        for(int i_formula = 0 ; i_formula < int(__variables_list__.size()) ; i_formula++ ){
          __observable_values__[__variables_list__[i_formula]] = __formula_list__[i_formula]->EvalInstance();
        }
        __outfile_ttree__->Fill();
      }
    }

    __outfile_tfile__->WriteObject(__outfile_ttree__, "selectedEvents");
    __outfile_tfile__->Close();
    cout << INFO << __outfile_tfile__->GetName() << " has been written." << endl;
    delete __outfile_tfile__;

    __current_tfile__->Close();
    delete __current_tfile__;

  }

  exit(0);

}

void fill_variables_expression(){

  cout << INFO << "Filling variables expression" << endl;

  // __cuts_string__ = "isRecoChargedLepton == 1 && flagCC0pi == 1 && isRecoProton == 1";
  __cuts_string__ = "isRecoChargedLepton == 1 && isRecoProton == 1 && isRecoPip == 0 && isRecoPim == 0 && Npi0 == 0 &&  && Nother==0";


  {
    __variables_list__.emplace_back("true_Enu");
    __expressions_list__.emplace_back("Enu_true");
  }

  {
    __variables_list__.emplace_back("true_mode");
    __expressions_list__.emplace_back("Mode");
  }

  {
    __variables_list__.emplace_back("true_mu_pT");
    string __lepton_observable_name__ = "pmu_4mom";
    // string __nucleon_observable_name__ = "hm_pprot_4mom";
    string px_sqr = "TMath::Power(" + __lepton_observable_name__ + ".Px(),2)";
    string py_sqr = "TMath::Power(" + __lepton_observable_name__ + ".Py(),2)";
    __expressions_list__.emplace_back("TMath::Sqrt(" + px_sqr + " + " + py_sqr + ")");
  }

  {
    __variables_list__.emplace_back("reco_mu_pT");
    string __lepton_observable_name__ = "pmu_reco_4mom";
    // string __nucleon_observable_name__ = "hm_pprot_4mom";
    string px_sqr = "TMath::Power(" + __lepton_observable_name__ + ".Px(),2)";
    string py_sqr = "TMath::Power(" + __lepton_observable_name__ + ".Py(),2)";
    __expressions_list__.emplace_back("TMath::Sqrt(" + px_sqr + " + " + py_sqr + ")");
  }

  {
    __variables_list__.emplace_back("true_prot_pT");
    // string __lepton_observable_name__ = "pmu_4mom";
    string __nucleon_observable_name__ = "hm_pprot_4mom";
    string px_sqr = "TMath::Power(" + __nucleon_observable_name__ + ".Px(),2)";
    string py_sqr = "TMath::Power(" + __nucleon_observable_name__ + ".Py(),2)";
    __expressions_list__.emplace_back("TMath::Sqrt(" + px_sqr + " + " + py_sqr + ")");
  }

  {
    __variables_list__.emplace_back("reco_prot_pT");
    // string __lepton_observable_name__ = "pmu_4mom";
    string __nucleon_observable_name__ = "hm_pprot_reco_4mom";
    string px_sqr = "TMath::Power(" + __nucleon_observable_name__ + ".Px(),2)";
    string py_sqr = "TMath::Power(" + __nucleon_observable_name__ + ".Py(),2)";
    __expressions_list__.emplace_back("TMath::Sqrt(" + px_sqr + " + " + py_sqr + ")");
  }

  {
    __variables_list__.emplace_back("true_delta_pT");
    string delta_pt_x = "(pmu_4mom.Px() + hm_pprot_4mom.Px())";
    string delta_pt_y = "(pmu_4mom.Py() + hm_pprot_4mom.Py())";
    __expressions_list__.emplace_back("TMath::Sqrt( TMath::Power(" + delta_pt_x + ",2)" + " + TMath::Power(" + delta_pt_y + ",2) )");
  }

  {
    __variables_list__.emplace_back("reco_delta_pT");
    string delta_pt_x = "(pmu_reco_4mom.Px() + hm_pprot_reco_4mom.Px())";
    string delta_pt_y = "(pmu_reco_4mom.Py() + hm_pprot_reco_4mom.Py())";
    __expressions_list__.emplace_back("TMath::Sqrt( TMath::Power(" + delta_pt_x + ",2)" + " + TMath::Power(" + delta_pt_y + ",2) )");
  }

  {
    __variables_list__.emplace_back("true_delta_alphaT");
    string __lepton_observable_name__ = "pmu_4mom";
    string __nucleon_observable_name__ = "hm_pprot_4mom";
    string delta_pt_x = "(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px())";
    string delta_pt_y = "(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py())";
    string norm_delta_pt = "TMath::Sqrt( TMath::Power(" + delta_pt_x + ",2)" + " + TMath::Power(" + delta_pt_y + ",2) )";
    string norm_pmu = "TMath::Sqrt( TMath::Power(" + __lepton_observable_name__ + ".Px(),2) + TMath::Power(" + __lepton_observable_name__ + ".Py(),2) )";
    string cos_delta_alpha = "-( " + __lepton_observable_name__ + ".Px() * " + delta_pt_x + " + " + __lepton_observable_name__ + ".Py() * " + delta_pt_y + ")/(" + norm_delta_pt + " * " + norm_pmu + ")";
    __expressions_list__.emplace_back("TMath::ACos(" + cos_delta_alpha + ")");
  }

  {
    __variables_list__.emplace_back("reco_delta_alphaT");
    string __lepton_observable_name__ = "pmu_reco_4mom";
    string __nucleon_observable_name__ = "hm_pprot_reco_4mom";
    string delta_pt_x = "(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px())";
    string delta_pt_y = "(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py())";
    string norm_delta_pt = "TMath::Sqrt( TMath::Power(" + delta_pt_x + ",2)" + " + TMath::Power(" + delta_pt_y + ",2) )";
    string norm_pmu = "TMath::Sqrt( TMath::Power(" + __lepton_observable_name__ + ".Px(),2) + TMath::Power(" + __lepton_observable_name__ + ".Py(),2) )";
    string cos_delta_alpha = "-( " + __lepton_observable_name__ + ".Px() * " + delta_pt_x + " + " + __lepton_observable_name__ + ".Py() * " + delta_pt_y + ")/(" + norm_delta_pt + " * " + norm_pmu + ")";
    __expressions_list__.emplace_back("TMath::ACos(" + cos_delta_alpha + ")");
  }

  {
    __variables_list__.emplace_back("true_Enu_QE");
    string __lepton_observable_name__ = "pmu_4mom";
    string __nucleon_observable_name__ = "hm_pprot_4mom";
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
    __expressions_list__.emplace_back(E_nu_QE);
  }

  {
    __variables_list__.emplace_back("reco_Enu_QE");
    string __lepton_observable_name__ = "pmu_reco_4mom";
    string __nucleon_observable_name__ = "hm_pprot_reco_4mom";
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
    __expressions_list__.emplace_back(E_nu_QE);
  }

  {
    __variables_list__.emplace_back("true_Enu_Calo");
    string __lepton_observable_name__ = "pmu_4mom";
    string __nucleon_observable_name__ = "hm_pprot_4mom";
    string proton_mass = "938.272";
    string E_nu_mu_p = "(";
    E_nu_mu_p += "" + __lepton_observable_name__ + ".Energy()";
    E_nu_mu_p += " + ";
    E_nu_mu_p += "0.5*( " + __nucleon_observable_name__ + ".Px()*" + __nucleon_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Py()*" + __nucleon_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Pz()*" + __nucleon_observable_name__ + ".Pz() )/" + proton_mass;
    E_nu_mu_p += ")";
    __expressions_list__.emplace_back(E_nu_mu_p);
  }

  {
    __variables_list__.emplace_back("reco_Enu_Calo");
    string __lepton_observable_name__ = "pmu_reco_4mom";
    string __nucleon_observable_name__ = "hm_pprot_reco_4mom";
    string proton_mass = "938.272";
    string E_nu_mu_p = "(";
    E_nu_mu_p += "" + __lepton_observable_name__ + ".Energy()";
    E_nu_mu_p += " + ";
    E_nu_mu_p += "0.5*( " + __nucleon_observable_name__ + ".Px()*" + __nucleon_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Py()*" + __nucleon_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Pz()*" + __nucleon_observable_name__ + ".Pz() )/" + proton_mass;
    E_nu_mu_p += ")";
    __expressions_list__.emplace_back(E_nu_mu_p);
  }

  {
    __variables_list__.emplace_back("true_fermi_momentum");
    string __lepton_observable_name__ = "pmu_4mom";
    string __nucleon_observable_name__ = "hm_pprot_4mom";
    string norm_pmu_L = __lepton_observable_name__ + ".Pz()";
    string E_mu = __lepton_observable_name__ + ".Energy()";
    string norm_pprot_L = __nucleon_observable_name__ + ".Pz()";
    string E_prot = __nucleon_observable_name__ + ".Energy()";
    string proton_mass = "938.272";
    string neutron_mass = "939.565";
    // EXAMPLE FOR CARBON
    string binding_energy = "92.16"; // Binding energy of the C12 (nuclear definition)
    string neutron_separation_energy = "20.3"; // Taking the 1p3/2 subshell -> most populated
    string Mass_A = "(6*" + proton_mass + " + 6*" + neutron_mass + " - " + binding_energy + ")";
    string Mass_A_minus_1 = "( " + Mass_A + " - " + neutron_mass + " + " + neutron_separation_energy + " )";
    string p_fermi_T = "TMath::Sqrt( ";
    p_fermi_T += "TMath::Power(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px(),2)";
    p_fermi_T += " + TMath::Power(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py(),2)";
    p_fermi_T += " )";
    string X_quantity = "(" + Mass_A + " + " + norm_pmu_L + " + " + norm_pprot_L + " - " + E_mu + " - " + E_prot + ")";
    string p_fermi_L = "( ";
    p_fermi_L += "0.5*" + X_quantity + " - ( TMath::Power(" + p_fermi_T + ",2) + TMath::Power(" + Mass_A_minus_1 + ",2) )/(2*" + X_quantity + ")";
    p_fermi_L += " )";
    string fermi_momentum = "TMath::Sqrt( TMath::Power(" + p_fermi_T + ",2) + TMath::Power( " + p_fermi_L + ",2) )";
    __expressions_list__.emplace_back(fermi_momentum);
  }

  {
    __variables_list__.emplace_back("reco_fermi_momentum");
    string __lepton_observable_name__ = "pmu_reco_4mom";
    string __nucleon_observable_name__ = "hm_pprot_reco_4mom";
    string norm_pmu_L = __lepton_observable_name__ + ".Pz()";
    string E_mu = __lepton_observable_name__ + ".Energy()";
    string norm_pprot_L = __nucleon_observable_name__ + ".Pz()";
    string E_prot = __nucleon_observable_name__ + ".Energy()";
    string proton_mass = "938.272";
    string neutron_mass = "939.565";
    // EXAMPLE FOR CARBON
    string binding_energy = "92.16"; // Binding energy of the C12 (nuclear definition)
    string neutron_separation_energy = "20.3"; // Taking the 1p3/2 subshell -> most populated
    string Mass_A = "(6*" + proton_mass + " + 6*" + neutron_mass + " - " + binding_energy + ")";
    string Mass_A_minus_1 = "( " + Mass_A + " - " + neutron_mass + " + " + neutron_separation_energy + " )";
    string p_fermi_T = "TMath::Sqrt( ";
    p_fermi_T += "TMath::Power(" + __lepton_observable_name__ + ".Px() + " + __nucleon_observable_name__ + ".Px(),2)";
    p_fermi_T += " + TMath::Power(" + __lepton_observable_name__ + ".Py() + " + __nucleon_observable_name__ + ".Py(),2)";
    p_fermi_T += " )";
    string X_quantity = "(" + Mass_A + " + " + norm_pmu_L + " + " + norm_pprot_L + " - " + E_mu + " - " + E_prot + ")";
    string p_fermi_L = "( ";
    p_fermi_L += "0.5*" + X_quantity + " - ( TMath::Power(" + p_fermi_T + ",2) + TMath::Power(" + Mass_A_minus_1 + ",2) )/(2*" + X_quantity + ")";
    p_fermi_L += " )";
    string fermi_momentum = "TMath::Sqrt( TMath::Power(" + p_fermi_T + ",2) + TMath::Power( " + p_fermi_L + ",2) )";
    __expressions_list__.emplace_back(fermi_momentum);
  }


}


void hook_tformula_on_input_tree(){

  for(auto const &formula : __formula_list__){
    delete formula;
  }
  __formula_list__.clear();

  for(int i_formula = 0 ; i_formula < int(__variables_list__.size()) ; i_formula++ ){
    __formula_list__.emplace_back( new TTreeFormula(
      __variables_list__[i_formula].c_str(),
      __expressions_list__[i_formula].c_str(),
      __current_ttree__)
    );
  }

}





////////////////////////////////////
//
// TOOLS
//
////////////////////////////////////
vector<string> get_list_of_files_in_directory(string folder_path_, string file_format_){

  vector<string> files_list;

  DIR* dirp = opendir(folder_path_.c_str());
  struct dirent * dp;
  while ((dp = readdir(dirp)) != NULL) {
    if( do_string_has_valid_format(dp->d_name, file_format_) ){
      files_list.emplace_back(dp->d_name);
    }
  }
  closedir(dirp);

  for(int i_file = 0 ; i_file < int(files_list.size()) ; i_file++){
    if(files_list[i_file] == "."
    or files_list[i_file] == ".."){
      files_list.erase(files_list.begin() + i_file);
      i_file--; // let the for loop do ++
    }
  }

  return files_list;

}

bool do_string_has_valid_format(string string_to_test_, string string_format_){

  vector<string> string_pieces = split_string(string_format_, "*");

  for(int i_piece = 0 ; i_piece < int(string_pieces.size()) ; i_piece++){

    if(string_pieces[i_piece].empty()) continue;

    std::vector<std::string> slices = split_string(string_to_test_, string_pieces[i_piece]);
    if(slices.size() == 1){
      return false;
    } else{

      vector<string> tail_slices;
      for(int i_slice = 0 ; i_slice < int(slices.size()) ; i_slice++){
        if(slices[i_slice].empty()) continue;
        tail_slices.emplace_back(slices[i_slice]);
      }
      string_to_test_ = join_vector_string(tail_slices, string_pieces[i_piece]);

    }
  }

  return true;

}


std::string join_vector_string(std::vector<std::string> string_list_, std::string delimiter_, int begin_index_, int end_index_) {

  std::string joined_string;
  if(end_index_ <= 0 and int(string_list_.size()) > fabs(end_index_)) end_index_ = int(string_list_.size()) + end_index_;

  for(int i_list = begin_index_ ; i_list < end_index_ ; i_list++){
    if(not joined_string.empty()) joined_string += delimiter_;
    joined_string += string_list_[i_list];
  }

  return joined_string;
}


std::vector<std::string> split_string(std::string input_string_, char delimiter_) {
  vector<string> output_vector;
  stringstream input_stringstream(input_string_); // Turn the string into a stream.
  string tokken;

  while(getline(input_stringstream, tokken, delimiter_)) {
    output_vector.push_back(tokken);
  }

  return output_vector;
}


std::vector<std::string> split_string(std::string input_string_, string delimiter_){

  vector<string> output_splited_string;

  const char *src = input_string_.c_str();
  const char *next = src;

  string out_string_piece = "";

  while ((next = strstr(src, delimiter_.c_str())) != NULL) {
    out_string_piece = "";
    while (src != next){
      out_string_piece += *src++;
    }
    output_splited_string.emplace_back(out_string_piece);
    /* Skip the delimiter_ */
    src += delimiter_.size();
  }

  /* Handle the last token */
  out_string_piece = "";
  while (*src != '\0')
  out_string_piece += *src++;

  output_splited_string.emplace_back(out_string_piece);

  return output_splited_string;

}


bool do_tfile_is_valid(TFile *input_tfile_, bool check_if_writable_){

  if(input_tfile_ == nullptr){
    if(verbosity_level >= 1) cerr << ERROR << "input_tfile_ is a nullptr" << endl;
    return false;
  }

  if(not input_tfile_->IsOpen()){
    if(verbosity_level >= 1) cerr << ERROR << "input_tfile_ = " << input_tfile_->GetName() << " is not opened." << endl;
    if(verbosity_level >= 1) cerr << ERROR << "input_tfile_->IsOpen() = " << input_tfile_->IsOpen() << endl;
    return false;
  }

  if(check_if_writable_ and not input_tfile_->IsWritable()){
    if(verbosity_level >= 1) cerr << ERROR << "input_tfile_ = " << input_tfile_->GetName() << " is not writable." << endl;
    if(verbosity_level >= 1) cerr << ERROR << "input_tfile_->IsWritable() = " << input_tfile_->IsWritable() << endl;
    return false;
  }

  return true;

}


void display_loading(int current_index_, int end_index_, string title_, bool force_display_) {

  int percent = int(round(double(current_index_) / end_index_ * 100.));
  if(force_display_ or current_index_ >= end_index_-1) {
    if(__last_displayed_value__ != -1) clog << "\r" << title_ << " : " << 100 << "%" << endl;
    __last_displayed_value__ = -1;
    return;
  }
  if(__last_displayed_value__ == -1 or __last_displayed_value__ < percent) {
    __last_displayed_value__ = percent;
    clog << "\r" << title_ << " : " << percent << "%" << flush << "\r";
  }

}
