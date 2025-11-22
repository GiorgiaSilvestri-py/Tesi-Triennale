import ROOT
import pandas as pd
import vector

def pdg_to_Name(pdgId):
    #map PDGs to particle names; 
    #return PDG value if not in dictionary
    pdgNames = { 
                11:"e-",  -11:"e+",    13:"mu-",  -13:"mu+",   15:"tau-",  -15:"tau+",
                12:"ve",  -12:"~ve",   14:"vmu",  -14:"~vmu",  16:"vtau",  -16:"~vtau",
                21:"g",    22:"γ",     23:"Z",     24:"W+",    -24:"W-",    25:"H",
                6:"t",    -6:"tbar",   5:"b",     -5:"bbar",   1:"d",      -1:"dbar",
                2:"u",    -2:"ubar",   3:"s",     -3:"sbar",   4:"c",      -4:"cbar"
               }
    return pdgNames.get(pdgId, f"PDG({pdgId})")

def build_4vec(pt, eta, phi, mass):
    return vector.obj(pt=pt, eta=eta, phi=phi, mass=mass)

def read_file(filename, maxEvents=50):
    file = ROOT.TFile.Open(filename)
    t = file.Get("Events")
    if not t:
        print("No Events tree in", filename)
        return pd.DataFrame()

    rows = []
    for i in range(min(t.GetEntries(), maxEvents)):

        print(f"----------- Event {i} -----------")

        t.GetEntry(i)

        e_vec = mu_vec = build_4vec(0,0,0,0)
        d1_vec = d2_vec = build_4vec(0,0,0,0)
        decay = "No Higgs"

        #LHEPart
        for j in range(int(t.nLHEPart)):
            pid = int(t.LHEPart_pdgId[j])       #get particle ID
            status = int(t.LHEPart_status[j])   #get particle status

            print(pid, status)

            if status == 1 :  # final state particle

                if abs(pid) == 11:
                    e_vec = build_4vec(t.LHEPart_pt[j], t.LHEPart_eta[j],
                                       t.LHEPart_phi[j], t.LHEPart_mass[j])

                    print(f"e vector = {e_vec}")

                elif abs(pid) == 13:
                    mu_vec = build_4vec(t.LHEPart_pt[j], t.LHEPart_eta[j],
                                        t.LHEPart_phi[j], t.LHEPart_mass[j])

                    print(f"mu vector = {mu_vec}")

        #Higgs decay - GenPart
        daughters, kin = [], []

        for j in range(int(t.nGenPart)):

            if int(t.GenPart_pdgId[j]) == 25:

                #get daughters' indices
                child_indices = [k for k in range(int(t.nGenPart))
                                 if int(t.GenPart_genPartIdxMother[k]) == j]

                #work on particles pair decays
                if len(child_indices) > 1:
                    for k in sorted(child_indices)[:2]:
                        daughters.append(pdg_to_Name(int(t.GenPart_pdgId[k])))
                        kin.append([t.GenPart_pt[k], t.GenPart_eta[k],
                                    t.GenPart_phi[k], t.GenPart_mass[k]])

                    print(f"daughters = {daughters} // kinematic info = {kin}")
                    break

        if len(kin) >= 2:
            decay = " + ".join(daughters)
            d1_vec = build_4vec(*kin[0])
            d2_vec = build_4vec(*kin[1])

        #add row to df
        rows.append({
            "event": i,
            "e_vec": e_vec,
            "mu_vec": mu_vec,
            "d1_vec": d1_vec,
            "d2_vec": d2_vec,
            "decay": decay
        })

    return pd.DataFrame(rows)

def expand_columns_vectors(df):
    vector_cols = []
    for col in df.columns:
        # controlla se la colonna contiene vector.obj
        if df[col].apply(lambda v: hasattr(v, "pt")).all():
            vector_cols.append(col)
            df[f"{col}_pt"]   = df[col].apply(lambda v: v.pt)
            df[f"{col}_eta"]  = df[col].apply(lambda v: v.eta)
            df[f"{col}_phi"]  = df[col].apply(lambda v: v.phi)
            df[f"{col}_mass"] = df[col].apply(lambda v: v.mass)

    df = df.drop(columns=vector_cols)

    return df


def main():
    df = read_file("nanoLatino_ZH_H_ZToLL_ZL__part0.root", maxEvents=50)
    df.to_pickle("events.pkl")

    print("Showing Pandas dataframe")
    print(df.head())

    rdf = expand_columns_vectors(df)
    print("Showing Expanded Pandas dataframe")
    print(rdf.head())

    #Pd dataframe --> ROOT df
    root_df = ROOT.RDF.FromPandas(rdf)
    root_df.Snapshot("Events", "file.root")

if __name__ == "__main__":
    main()