import root_pandas as rp
import numpy as np
import matplotlib
from matplotlib import gridspec
import matplotlib.pyplot as plt
import pandas as pd


print( "loading data frame named df.")
df = rp.read_root("data/pan_process_Sept.root")
print("finished loading data frame")


print("unc_mc_process and scales are dictionaries")
# uncertainty on process
unc_mc_process = { "WW": .05, "DY": .03, "TT": 0.05, "ZZ": 0.1, "WZ":.1 }
# process ratio scalings

scales = {"WW":       19.7e3 * 59.8 /  10000431.0 ,\
          "DY":       19.7e3 * 3531.9 / 30459500. ,\
          "TT":       19.7e3 * 25.81 / 12011428.,\
          "ZZ":       19.7e3 * 9.03 / 9799908.,\
          "WZ_2l2q":  19.7e3 * 2.24 / 3215990,\
          "WZ_3ln":  19.7e3 * 1.07 / 2017979. ,\
          "WZ":       19.7e3 * 1.07 / 2017979., "WJ": 1, "Da": 1 }


colors = {"WW": (0.993248, 0.906157, 0.143936),\
          "DY": (0.344074, 0.780029, 0.397381),\
          "TT": (0.121380, 0.629492, 0.531973),\
          "ZZ": (1. ,  0.85490196,  0.7254902 ),\
          "WZ_2l2q": (0.165117, 0.467423, 0.558141),\
          "WZ_3ln": (0.165117, 0.467423, 0.558141),\
          "WZ": (0.165117, 0.467423, 0.558141),\
          "WJ": (0.267004, 0.004874, 0.329415)}


scale_data = 5.27/19.7
print( "scale_data = (5.27/19.7)")

def create_table( data, round_digit ):
    for flavor in data.keys():
        print flavor
        print "\t",{ process : round(data[flavor][process], round_digit) for process in data[flavor].keys() }

def combine_unc( data ):
    comb_unc = {}
    for process in data[data.keys()[0]].keys():
        comb_unc[process] = 0
    for flavor in data.keys():
        for process in data[flavor].keys():
            comb_unc[process] += pow( data[flavor][process], 2 )
    print { process : round( pow( comb_unc[process], .5) , 2) for process in comb_unc.keys()}



def cuts_ana( df ):
  """
    Basic cuts analysis
    cuts_ana( df ) 
    df: data frame
  """
  initial_cuts = (df.mll > 12 ) & (df.qT > 30) & (df.METProj > 20) & (df.numbExtraLep == 0)
  same_flavor_cut = (df.lep_Type < 0)
  dif_flavor_cut = (df.lep_Type > 0)
  z_peak_cut = (df.mll > 106) | (df.mll < 76)
  met_proj_cut = (df.METProj > 45)

  basic_sf_0j_cuts = initial_cuts & (df.qT > 45) &same_flavor_cut & met_proj_cut & z_peak_cut& (df.numb_jets <= 1)
  basic_df_0j_cuts = initial_cuts & dif_flavor_cut &(df.numb_jets <= 1)
  return pd.concat( [df[basic_sf_0j_cuts], df[basic_df_0j_cuts]] )


def WW_ana( df ): 
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets ==  0) & ( df.metMod > 50 )
  dif_flavor_cut = (df.lep_Type > 0)

  basic_df_0j_cuts = initial_cuts & dif_flavor_cut 
  return df[basic_df_0j_cuts]


def TT_ana( df ): 
  initial_cuts = (df.mll > 50 ) & (df.qT > 30) & (df.metMod > 60) & (df.numbExtraLep == 0) & (df.numb_jets > 0)
  dif_flavor_cut = (df.lep_Type > 0)
  
  basic_df_0j_cuts = initial_cuts & dif_flavor_cut 
  return df[basic_df_0j_cuts] 

def DY_ana( df ): 
  initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets == 0) & (df.metMod < 80)

  return df[ initial_cuts ]

def Z_tt_ana( df ): 
    initial_cuts = (df.mll > 50 )  & (df.numbExtraLep == 0) & (df.numb_jets ==  0) & ( df.qT < 10 ) & (df.lep1_pt > 30) & (df.lep2_pt >20)
    
    dif_flavor_cut = (df.lep_Type > 0)
    
    m_1 = (df.mll < 120) & (df.mll > 60)
    met_1 = (df.metMod > 30) & (df.metMod < 50)

    return df[ initial_cuts & met_1 & m_1 & dif_flavor_cut ] 
###############################################
#### Create kinematic histograms
def create_kinematic_hist(df_mc, df_data, prefix=""):
  """
  Creates all the basic histograms you'll ever need:
  create_kinematic_hist(df)
  """
  plots = []
  ratio_plots = []

  range = (0,250)
  bins  = 100

  features = ['sub_lep_eta', 'nBJet', 'HT',
       'numb_jets', 'lead_lep_eta', 'lep1_pt',
       'jetPt1', 'lep2_pt',
       'jetPt3', 'jetPt2', 'metMod', 'dPhiLL', 'METProj',
       'qT', 'dPhiLLJet', 'MET_phi',
       'dPhiLLMET', 'lep_Type', 'METProj_parr', 'lep3_pt',
       'mll']

  for feature in features:
    if ('Phi' in feature) or ('phi' in feature) or ('eta' in feature):
      range = (0,4)
      bins = 100
    elif ('numb_jet' in feature) or ('nBJet' in feature):
      range = (-.5, 9.5)
      bins = 10
    else:
      range = (0,250)
      bins  = 100

    print feature
    bins_mc = bin_df( df_mc, feature, range=range, scales=scales_sept) 
    bins_data = bin_df( df_data, feature, range=range, scales={'Da':1}, lumi_scale=1) 
    figs, ax = full_plot(bins_mc, bins_data, )#logy=False, y_range=(0,200))
    figs.savefig('plots/'+prefix+"_"+feature+".png")

  return plots
##############################################
###   Below is Depricated 
'''
    in the future use plt.histo( [a,b], stacked=True, weights=[ a.shape[0]*[1], b.shape[0]*[.6] ])

    **There is no need for the bin_df step and rect placement.**
'''

###############################################
#### Histogramming and Binning
def bin_df( df, binned_feature, scales=scales, range=None, bins=100, lumi_scale=1):
  """
  Binning stuffs:
  bin_df( df, binned_feature, range=None,bins=100, lumi_scale=1)
  """

  bins_rf = {}
 
  new_df = 0 
  weights = None
  if "weights" in df.keys():
    weights = df.weights.values*lumi_scale
  for process in scales.keys():
      if process not in df.keys():
          new_df += 1
          continue
      bins_rf[process] = list(np.histogram(df[ (df[process]==1) ][binned_feature] , bins=bins, range=range, weights=weights*scales[process]))
      bins_rf[process][0] = bins_rf[process][0] * scales[process]
      bins_rf[process].append( bins_rf[process][1][0:-1] + ( bins_rf[process][1][1:] - bins_rf[process][1][0:-1] ) / 2. )
      bins_rf[process].append( np.histogram(df[ (df[process]==1) ][binned_feature] , bins=bins, range=range, weights=(weights**2)*scales[process])[0] )

  if new_df != 0:
      #print "new version"
      for process in scales.keys():
        bins_rf[process] = list(np.histogram(df[ (df["process"]==process) ][binned_feature] , bins=bins, range=range,  weights=weights*scales[process]))
        bins_rf[process][0] = bins_rf[process][0] * scales[process]*lumi_scale
        bins_rf[process].append( bins_rf[process][1][0:-1] + ( bins_rf[process][1][1:] - bins_rf[process][1][0:-1] ) / 2. ) 
        bins_rf[process].append( np.histogram(df[ (df[process]==1) ][binned_feature] , bins=bins, range=range, weights=(weights**2)*scales[process])[0] )

  return bins_rf




def plot_hist( bins, processes=[ "WW", "TT", "WZ", "ZZ", "DY"], x_range=None, y_range=None, title=None, y_label=None, color=colors, logy=True, x_minor_ticks=True, ax=None):
  """
  Histogramming stuffs
  plot_hist( bins, processes=[ "WW", "TT", "WZ", "ZZ", "DY"], x_range=None, y_range=None, title=None, y_label=None, color=colors, logy=True, x_minor_ticks=True)
  """
  if ax == None: fig, ax = plt.subplots(figsize=(11, 5))

  rect = []
  sum_bins = np.zeros(len(bins[processes[0]][0]))

  previous_process = []


  for process in processes:
#########
    for key in color.keys():
      if process in key:
        process_ = key
        if process_ not in previous_process:
            previous_process.append( process_ )
        else:
            continue         
########
        bottom = sum_bins
        #print process_
        rect.append(ax.bar( bins[process_][1][:-1], bins[process_][0],
                      bins[process_][1][1]-bins[process_][1][0] , color = colors[process_],
                      edgecolor = colors[process_], bottom=bottom ))
        sum_bins +=bins[process_][0]

#  for process in processes:
#########
#    for key in color.keys():
#      if process in key:
#        process_ = key
#        previous_process.append( process_ )
########
#    bottom = sum_bins
#    rect.append(ax.bar( bins[process][1][:-1], bins[process][0],
#                      bins[process][1][1]-bins[process][1][0] , color = colors[process],
#                      edgecolor = colors[process], bottom=bottom ))
#    sum_bins +=bins[process][0] 
#    rect.append(ax.bar( bins[process_][1][:-1], bins[process_][0],
#                      bins[process_][1][1]-bins[process_][1][0] , color = colors[process_],
#                      edgecolor = colors[process_], bottom=bottom ))
#    sum_bins +=bins[process_][0] 

  #Configurables
  if logy == True: ax.set_yscale("log", nonposy='clip')
  if x_range!=None: ax.set_xlim(x_range)
  if y_range==None: ax.set_ylim( bottom=1,  top= sum_bins.max()*30.)
  elif type(y_range)==tuple: ax.set_ylim( y_range )

  if y_label != None:
      plt.ylabel(y_label, position=(-0.25, 1.), va='top', ha='right', fontsize=18, fontname='Arial')
  if title != None:
      plt.xlabel( title, position=(1., 0.), va='bottom', ha='right', fontsize=18, fontname='Arial')
 

  ####################################
  #Add minor tick marks to the x-axis
#  if x_minor_ticks == False:
#      loc = matplotlib.ticker.MultipleLocator(base=1) # this locator puts ticks at regular intervals
#      ax.xaxis.set_major_locator(loc)
#  else:
#      loc = matplotlib.ticker.MultipleLocator(base=(bins[process][1][-1] - bins[process][1][0])/25.) # this locator puts ticks at regular intervals
#      ax.xaxis.set_minor_locator(loc)
  
###################################
  ax.yaxis.set_tick_params(length=10)
  ax.yaxis.set_tick_params(which='minor',length=5)
  ax.xaxis.set_tick_params(length=10)
  ax.xaxis.set_tick_params(which='minor',length=5)
      
      
  ax.xaxis.labelpad = 20
  ax.yaxis.labelpad = 20
  
  plt.tight_layout()

  return ax, rect, processes


   
def plot_errbar( bins, process='Da', ax=None ):
  """
  Plot errbar
  plot_errbar( bins, process='Da' )
  bins: dictionary of list containing bin contents, edges and middles
  """
  if ax == None: plotline, caplines, barlinecols = plt.errorbar( bins[process][2], bins[process][0], yerr=np.sqrt(bins[process][0]), ecolor='black',color="black",fmt="o", label="data" )
  else: plotline, caplines, barlinecols = ax.errorbar( bins[process][2], bins[process][0], yerr=np.sqrt(bins[process][0]), ecolor='black',color="black",fmt="o", label="data" )

  return plotline, caplines, barlinecols 


def plot_ratio( bins_1, bins_2, y_label=None, x_label=None, ax=None):
  """
  Plot ratio plots
  plot_ratio( bins_1, bins_2, y_label=None, x_label=None, ax=None):
  bins_(1/2): list of three numpy arrays [ bins contents, bin edges, bin centers ]
  """
  tot_1 = np.zeros( bins_1[ bins_1.keys()[0] ][0].shape[0] )
  for process in bins_1.keys():
    tot_1 += bins_1[process][0]

  key = ''
  tot_2 = np.zeros( bins_2[ bins_2.keys()[0] ][0].shape[0] )
  for process in bins_2.keys():
    tot_2 += bins_2[process][0]
    key = process

  if ax==None : plotline, caplines, barlinecols =  plt.errorbar( bins_2[key][2], tot_1 / tot_2, ecolor='black',color="black",fmt="o" )
  else: 
    plotline, caplines, barlinecols =  ax.errorbar( bins_2[key][2], tot_1 / tot_2, ecolor='black',color="black",fmt="o" )
    ax.set_ylim( bottom=0.5,  top= 1.5)
    ax.set_ylabel('mc/data')

  return plotline, caplines, barlinecols




def full_plot(bins_mc, bins_data,  processes=[ "WW", "TT", "WZ", "ZZ", "DY"], x_range=None, y_range=None, title=None, y_label=None, color=colors, logy=True, x_minor_ticks=True):

  fig, ax = plt.subplots(2, figsize=(11,8), sharex=True, gridspec_kw = {'height_ratios':[3, 1]})

  hist_stuff = plot_hist( bins_mc, processes=[ "WW", "TT", "WZ", "ZZ", "DY"], x_range=x_range, y_range=y_range, title=title, y_label=y_label, color=color, logy=logy, x_minor_ticks=x_minor_ticks, ax=ax[0])
  err_stuff  = plot_errbar( bins_data, ax=ax[0])
  ax[0].legend(hist_stuff[1] + [err_stuff], hist_stuff[2] + ['Da'], frameon=False, fontsize="x-large", numpoints=1)

  plot_ratio( bins_mc, bins_data, ax=ax[1])

  fig.subplots_adjust(hspace=0.1)

  return fig, ax
