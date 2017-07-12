typedef ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<float> > LV;
typedef std::vector<LV> LVs;

// The magic line.
// std::vector<ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<float> >>* lep_p4 = (std::vector<ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<float> >>*) t->GetBranch("lep_p4")->GetLeaf("lep_p4_")->GetValuePointer()
LVs* getLVPointers(TTree* t, TString branch_name, int entry=-1)
{
    if (entry >= 0)
        t->GetEntry(entry);
    // The following returns (void *)
    return (LVs*) t->GetBranch(branch_name)->GetLeaf(branch_name + "_")->GetValuePointer();
}

// The magic line.
// std::vector<ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<float> >>* lep_p4 = (std::vector<ROOT::Math::LorentzVector< ROOT::Math::PxPyPzE4D<float> >>*) t->GetBranch("lep_p4")->GetLeaf("lep_p4_")->GetValuePointer()
const LVs& getLVs(TTree* t, TString branch_name, int entry=-1)
{
    if (entry >= 0)
        t->GetEntry(entry);
    // The following returns (void *)
    return (*((LVs*) t->GetBranch(branch_name)->GetLeaf(branch_name + "_")->GetValuePointer()));
}

std::vector<bool> getBoolVec(TTree* t, TString branch_name, int entry=-1)
{
    if (entry >= 0)
        t->GetEntry(entry);
    // The following returns (void *)
    return (*((std::vector<bool>*) t->GetBranch(branch_name)->GetLeaf(branch_name)->GetValuePointer()));
}

std::vector<int> getIntVec(TTree* t, TString branch_name, int entry=-1)
{
    if (entry >= 0)
        t->GetEntry(entry);
    // The following returns (void *)
    return (*((std::vector<int>*) t->GetBranch(branch_name)->GetLeaf(branch_name)->GetValuePointer()));
}

int d_entry=11254;