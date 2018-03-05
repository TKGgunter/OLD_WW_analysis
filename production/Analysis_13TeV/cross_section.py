from prep_ana_II import *

scales = scales
lumi_amount = lumi_amount
def Xs_ingredients(df, da, scales=scales, yield_style=None):
    yield_df = process_yields(df, df_da, scales=scales, yield_style=None)

    Br = (3 * 0.108)**2.
    acceptance = 0.18
    lumi = float(lumi_amount) 

    background = 0
    signal = 0
    total_sig = 0
    data = 0

    return {"Br": Br, "Acc": acceptance, "Lumi": lumi, "Data": data,\
            "Bkg": background, "Sig": signal, "Tot_Sig": total_sig}
    

def calc_Xs(df, da, scales=scales, yield_style=None, fiducial=False, from_mc=False, **kwgs):
    if type(kwgs) == None:
        ingredients = Xs_ingredients(df, da, scales, yield_style)
    else:
        ingredients = kwgs

    Br = ingredients["Br"]
    Acc = ingredients["Acc"]
    Lumi = ingredients["Lumi"]
    Data = ingredients["Data"]
    Bkg = ingredients["Bkg"]
    Sig = ingredients["Sig"]
    Tot_sig = ingredients["Tot_Sig"]

    Xs = 0
    if fiducial == True:
        if from_mc == True:
            Xs = Sig / (Lumi * Sig / Tot_sig)
        if from_mc == False:
            Xs = (Data - Bkg) / (Lumi * Sig / Tot_sig)
    if fiducial == False:
        if from_mc == True:
            Xs = Sig / (Lumi * Acc * Br * Sig / Tot_sig)
        if from_mc == False:
            Xs = (Data - Bkg) / (Lumi * Acc * Br * Sig / Tot_sig)
   return Xs 
