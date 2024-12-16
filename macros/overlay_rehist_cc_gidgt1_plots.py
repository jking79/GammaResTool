#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *
#import cmsstyle as CMS

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t ):

    first = True
    #dofit = True
    dofit = False
    #sxtal = True
    sxtal = False
    paramn = []
    parnerror = []
    paramc = []
    parcerror = []
    params = []
    parserror = []
    thebinmid = []
    thebinerror = []
    f1 = []
    h1 = []
    n = 0


    setTDRStyle()
    gROOT.SetBatch(True)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)
    gStyle.SetPadTopMargin(0.09)
    gStyle.SetLabelOffset(0.005, "XYZ")
    #c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
    c1 = TCanvas( 'c1', 'canvas' , 300, 50, 1200, 900 )
    c1.cd()
    #if layout['logx'] : c1.SetLogx()
    if layout['logx'] : gStyle.SetOptLogx(1)
    if layout['logy'] : gStyle.SetOptLogy(1)
    if layout['logy'] : c1.SetLogy()
    c1.SetGridx(1)
    c1.SetGridy(1)
    c1.Update()
    c1.Draw()

    legend = TLegend(l[0],l[1],l[2],l[3]);
#    legend.SetHeader(layout['legtitle'], '')
    legend.SetName('legend')
    gStyle.SetLegendFont(42)
    gStyle.SetLegendTextSize(0.03)
    gStyle.SetLabelSize(0.03, "XYZ")
    gStyle.SetTitleSize(0.04, "XYZ")
    gStyle.SetTitleXOffset(1.3)
    gStyle.SetTitleYOffset(1.8)
    legend.SetBorderSize(0)
    #legend.SetLineColor(kBlack)

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    if dofit :
        ns=str(n)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',75,2250,2)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',75,500,2)
        #hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] ) )',100,2250,2)
        #hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] ) )',100,700,2)
        hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',75,3000,3)
        #hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',75,375,3)
        #hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',100,750,0,3)
        #hfit = TF1('hfits','sqrt((([0]*[0])/(x*x))+(2*[1]*[1]))',6,100,2)
        hfit.SetParName(0,'N')
        hfit.SetParameter(0,40.0)
        #hfit.SetParLimits(0,0,50)
        hfit.SetParLimits(0,0.0,100.0)
        #hfit.SetParameter(0,5)
        #hfit.SetParLimits(0,0,10)
        hfit.SetParName(1,'C')
        hfit.SetParameter(1,0.1)
        hfit.SetParLimits(1,0.0,1.0)
        #hfit.SetParLimits(1,0.01,10.0)
        #hfit.SetParLimits(1,0.02,1.0)
        #hfit.SetParameter(1,0.05)
        #hfit.SetParLimits(1,0.001,1.0)
        hfit.SetParName(2,'S')
        hfit.SetParameter(2,5.0)
        hfit.SetParLimits(2,0.0,25.0)

    mg = TMultiGraph();

    #thebinmid.append([  88.01, 111.9, 136.4, 161.4, 196.6, 247.8, 298.5, 348.9, 419.7, 528.9, 656.9, 801.5,   990, 1487.5,  1975])
    #thebinerror.append([ 7.09, 7.156, 7.154, 7.172, 14.28, 14.39, 14.43,  14.4, 28.53, 35.37, 40.58, 43.16, 38.06,      0,     0])
    #thebinmid.append([  90.43, 112.6, 137.1, 161.8, 196.8, 247.1, 297.6, 347.9, 417.6, 528.2,   663, 825.1,  1022,   1320,  1975])
    #thebinerror.append([6.521, 7.156, 7.192, 7.192, 14.22, 14.29, 14.33, 14.34,  28.4, 35.39, 42.34, 53.72, 63.78,  51.57,     0])
    #thebinmid.append([  94.04,   115, 137.7, 162.2, 198.1, 247.7, 297.7, 347.9,   418, 528.4, 664.2, 833.9,  1064,   1387,  1779])
    #thebinerror.append([5.056, 6.691, 7.158, 7.201, 14.32, 14.31, 14.33, 14.36, 28.43, 35.41, 42.68, 56.53,  84.3,  81.16, 65.90])

    for histname, tree, infile, lego  in hist_list :
    
        f1.append(TFile.Open(infile))
        if tree == '' : hist = histname
        else : hist = tree+'/'+histname
 
        orighist = f1[n].Get(hist)
        #nbins = 15 #orighist.GetNbinsX()
        #low = 0.0 #orighist.GetMinimumBin()
        #high = 100.0 #orighist.GetMaximumBin()
        #bscale = (high-low)/(10*nbins)
        #print('Find Bins',nbins,low,high,bscale) 
        htitle = 'hist' + str(n)

        #lenmybins = len(mybins) - numExcludeBins
        lenmybins = int(orighist.GetNbinsX())
        #lenmybins = 16
        #myhist = TH1D(htitle,"",150,low,high)
        #h1.append(TGraphErrors(lenmybins))
        h1.append(TGraphErrors(lenmybins))
        #h1.append(TH2D(htitle,"",(10*nbins),low,high,240,0,2.4))
        #h1.append(TH1F(htitle,"",len(mybins)-1,mybins))

        for bn in range( 0,lenmybins+1):
            #binval = float(orighist.GetBinContent(bn))*binwidths[bn-1]*2
            binval = float(orighist.GetBinContent(bn))
            #binerr = 0.0
            binerr = float(orighist.GetBinError(bn))
            binmid = float(orighist.GetBinCenter(bn))
            binwidth = float(orighist.GetBinWidth(bn))
            binstart = float(orighist.GetBinLowEdge(bn))
            #binmid = thebinmid[n][bn-1]  #float(orighist.GetBinCenter(bn))
            if binval > 0.0 :     
            #goodpoint = binerr < (binval/5.0) and binerr != 0
            #if goodpoint :
                if( sxtal ):  
                    print( '--- adjusting bin value' )
                    binval = binval/sqrt(2)
                h1[n].SetPoint(bn,binmid,binval)
                #ebin = int(binmid/bscale)+1
                widtherr = (binwidth/2)/sqrt(3)
                #widtherr = 0.0004
                #widtherr = (0.25)/sqrt(3)
                #widtherr = thebinerror[n][bn-1] #binwidths[bn-1]/sqrt(3)
                h1[n].SetPointError(bn,widtherr,binerr)
                print('Fill bin',bn,'at',binmid,'with',binval,'err width',widtherr,' err bin',binerr,'for:',binstart,'to',binstart+binwidth,'width',binwidth) 
                #print('Fill bin',bn,'at',binmid,'with',binval,'for',ebin,'with',binerr)  
        #h1[n].SetPoint(lenmybins,10000,0.001)
#       h1.append(f1[n].Get(hist))

#       print( f1 )
#       print( h1 )

        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #h1[n].SetMarkerStyle(6)
        #h1[n].SetTitle(layout['title'])
        #h1[n].GetXaxis().CenterTitle(True)
        #h1[n].GetXaxis().SetTitle(layout['xtitle'])
        #h1[n].GetYaxis().CenterTitle(True)
        #h1[n].GetYaxis().SetTitle(layout['ytitle'])
#       k = [kMagenta+2,kBlue+1,kAzure+4,kBlack,kYellow+1,kViolet+7,kOrange+7,kRed+2,kGreen+3, kGray]
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        #k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        k = [kGreen+2,kBlue+2,kRed+2,kAzure+4,kGreen+3,kOrange+7,kViolet+7,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+2,kGray+3,kGray+4,kBlack]
        #k = [kGreen+2,kBlue+2,kMagenta+2,kRed+2]
        h1[n].SetLineColor(k[n])
        #h1[n].SetLineWidth(1)
        #h1[n].SetLineStyle(9-n)
        if dofit : 
            hfit.SetLineColor(k[n])
            hfit.SetLineWidth(2)
            #hfit.SetLineStyle(9-n)
        #if dofit : hfit.SetLineColor(kRed+2)
        h1[n].SetMarkerColor(k[n])
        #msz = 0.2
        #msz = 0.8 # my standard
        msz = 1.6 # prsentaion plots
        #msz = 1.0
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        elif( n == 4 ) : h1[n].SetMarkerSize(msz+0.5)
        elif( n == 5 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)
        #h1[n].SetMinimum(y[1])
        #h1[n].SetMaximum(y[0])
#       h1[n].GetXaxis().SetRangeUser(200.0,1100.0)
        #h1[n].GetXaxis().SetRangeUser(x[0],x[1])
        #if layout['logx'] : h1[n].GetXaxis().SetMoreLogLabels()
        #if layout['logy'] : h1[n].GetYaxis().SetMoreLogLabels()
        #c1.cd()
        #c1.Update()
        if( first ) :
                #h1[n].Draw('ep')
                #h1[n].Draw('AP')
                if dofit : h1[n].Fit(hfit,'RE')
                #if dofit : h1[n].Fit(hfit,'REQ')
                first = False
        else :
                #h1[n].Draw('epSAME')
                #h1[n].Draw('SAME AP')
                if dofit : h1[n].Fit(hfit,'RE+')
                #if dofit : h1[n].Fit(hfit,'REQ+') 
        #c1.Update()
        #c1.Draw()

        if dofit : 
                 paramn.append(str(abs(hfit.GetParameter(0))))
                 paramc.append(str(abs(hfit.GetParameter(1))))
                 params.append(str(abs(hfit.GetParameter(2))))
                 pne = hfit.GetParError(0)
                 pce = hfit.GetParError(1)
                 pse = hfit.GetParError(2)
                 #print('Fit info',paramn[n],pne,paramc[n],pce)
                 print('Fit info',paramn[n],pne,paramc[n],pce,params[n],pse)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror.append(str(pne))
                 parcerror.append(str(pce))
                 parserror.append(str(pse))
                 #print( 'C: ' + param + ' +/- ' + parerror )
                 #lat_param = '#color['+str(k[n])+']{N : ' + paramn[0:4] + ' #pm ' + parnerror[0:4] + ' [ns]  C : ' + paramc[0:6] + ' #pm ' + parcerror[0:6] + ' [ns]}'
                 #lat.SetTextSize(0.03);
                 #lat.DrawLatex(t[3],t[4]-n*.035,lat_param);    
                 #c1.Modified()
                 #c1.Update()       
 
        legend.AddEntry(h1[n],lego,'epl');
        #c1.Update()
        #c1.Draw()
        n += 1

        #End of loop
    
    for h in range(0,n):
        mg.Add(h1[h])

    mg.Draw('AP')

    #mg.UseCurrentStyle()
    #mg.SetMarkerStyle(n+25)
    #mg.SetMarkerStyle(6)
    mg.SetTitle(layout['title'])
#+';X axis '+layout['xtitle']+';Y Axis '+layout['ytitle']+';')
    mg.GetXaxis().CenterTitle(True)
    mg.GetXaxis().SetTitle(layout['xtitle'])
    mg.GetYaxis().CenterTitle(True)
    mg.GetYaxis().SetTitle(layout['ytitle'])
    mg.SetMinimum(y[1])
    mg.SetMaximum(y[0])
#   mg.GetXaxis().SetRangeUser(200.0,1100.0)
    mg.GetXaxis().SetRangeUser(x[0],x[1])
    if layout['logx'] : mg.GetXaxis().SetMoreLogLabels()
    if layout['logy'] : mg.GetYaxis().SetMoreLogLabels()

    #mg.Draw('AP')
    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch

    gPad.Modified()
    #c1.SetGridx(1)
    #c1.SetGridy(1)
    #c1.cd()
    #c1.Update()
    #mg.Draw('AP')

    lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    #lat_cms = '#bf{CMS} #it{Work in Progress}' + ptitle[0]
    #lat_title = 'Run2018D 3206730-320824' #   7 fb^{-1} (#sqrt{s} = 13 TeV)'
    #lat_title = 'Run2018D 1Tier miniAOD'
    lat_title = ptitle[1]+' (13.6 TeV)'
    #lat_title = ptitle[1]+' (13.8 TeV)'
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (N/Eeff)^{2} + 2C^{2}'
    #lat.SetTextSize(0.045);
    #lat.SetTextFont(132);
    #lat.DrawLatex(0.48+max(0.,0.32-lat.GetXsize()),0.947,lat_string);
    lat.SetTextSize(0.045);
    lat.SetTextFont(42);
    lat.DrawLatex(0.16,0.93,lat_cms);
    #lat.DrawLatex(0.12,0.9325,lat_cms);
    #lat.DrawLatex(0.15,0.9425,lat_cms);
    #lat.DrawLatex(0.25,0.9325,lat_cms);
    #lat.SetTextSize(0.04);
    #lat.SetTextFont(42);
    #lat.DrawLatex(0.58,0.9325,lat_title);
    lat.DrawLatex((0.828-t[2]),0.93,lat_title);
    lat.SetTextSize(0.045);
    lat.DrawLatex(t[0],t[1],ptitle[2]);
    if dofit : 
        lat.SetTextSize(0.03);
        lat.DrawLatex(t[3],t[4]+.08,lat_form);
        #lat.SetTextAlign(12)
        #lat.DrawLatex(.74,0.0325,layout['xtitle']);
        for l in range(0,n):
            #lat_param ='#color['+str(k[l])+']{N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:4]+' [ns]  C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' }'
            lat_param =	'#color['+str(k[l])+']{'
            lat_param = lat_param + 'N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:3]+' [ns]   '
            lat_param = lat_param + 'S : '+params[l][0:3]+' #pm '+parserror[l][0:3]+' [ns]   '
            lat_param = lat_param + 'C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' [ns]}'
            lat.SetTextSize(0.03);
            lat.DrawLatex(t[3],t[4]-l*.035,lat_param);

    
    if layout['logx'] : c1.SetLogx()
    ####c1.BuildLegend()
    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.SaveAs( outname + '_' + date + '.C' )
    #c1.Show()
    c1.Close()

#from overlay_hist_defs_v2 import *
#from overlay_hist_defs_v3 import *

#dispho_2t_eg2018B_preplot_sp

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

#legtitle = 'Full dR '
#legtitle = 'dR < 2.5 '
#legtitle = 'dR > 2.5 '
legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run3' # (13TeV)'

gi_legtitle = 'Global Inclusive'
li_legtitle = 'Local Inclusive'
ls_legtitle = 'Local Same'
Ic_legtitle = ''
loc2_legtitle = ' LS'

xtitle = 'Energy [GeV]'
#xtitle = 'A_{eff}/#sigma_{n}'
#xtitle = '#Delta_{Run}'
#xtitle = 'E_{eff} [GeV]'
#xtitle = 'RecHit E [GeV]'
#xtitle = 'pCaloHit E_{eff} [GeV]'
#xtitle = 'pCaloHit energy [GeV]'
#xtitle = 'GeV'
#xtitle = '[ns]'
#xtitle = 'Max(R_{oot})'
#xtitle = 'Max(R_{oot}^{after})'
#xtitle = 'Max(R_{oot}^{before})'
#ytitle = '#sigma(t_{1}-t_{2}) [ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(t_{gen}-t_{reco}) [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
ytitle = ''
#htitle = 'A_{eff}/#sigma_{n} vs #sigma_{#delta_{t}}'
htitle = ''
#islogx = True
islogx = False
islogy = True
#islogy = False

gll_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
loc_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle + li_legtitle }
#s2legtitle = '#splitline{'+legtitle+loc2_legtitle+'}{Cali w/ TOF}'
#s2legtitle = '#splitline{'+legtitle+loc2_legtitle+'}{Cali w/out TOF}'
#s2legtitle = '#splitline{'+legtitle+loc2_legtitle+'}{Cali TOF Comp}'
s2_li_legtitle = '#splitline{'+legtitle+'}{'+li_legtitle+'}'
li_legtitle = legtitle+' '+li_legtitle
s2_gi_legtitle = '#splitline{'+legtitle+'}{'+gi_legtitle+'}'
gi_legtitle = legtitle+' '+gi_legtitle
Ic_legtitle = legtitle+' '+Ic_legtitle

loc2_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : li_legtitle }
glo2_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : gi_legtitle }
Ic_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : Ic_legtitle }
glo_layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle + gi_legtitle }
#---------------------------------------------------------------

layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }

#
##egres_r3_24_EG01_v12 = '../../res_files/egmres_24CD_Full_mixed_resplots/egres_Run2024CD_mixed_eg01_Full_140_v10_resplots.root'
##egres_r3_24_EG01_v12 = '../../res_files/egmres_24CD_Full_mixed_cali_resplots/egres_Run2024CD_mixed_eg01_Full_cali_140_v10_resplots.root'
#egres_r3_24_EG01_v12 = '../../res_files/egmres_24E_test_resplots/egres_Run2024E_test_140_v10_resplots.root'
##egres_r3_23d_part_EG1_v12 = '../../res_files/egmres_23D_partial_nojson_resplots/egres_Run2023D_SF_140_v10_resplots.root'
#egres_r3_23d_part_EG1_v12 = '../../res_files/egmres_23D_partial_nojson_cali_resplots/egres_Run2023D_SF_140_v10_resplots.root'

histfile12 = "ecal_cc_diag/egammares_diag_jetmet1_aod_tevjets_v3_gidgt1_rt1cc1_p12_kWDW_v21_diag.root"
histfile3 = "ecal_cc_diag/egammares_diag_jetmet1_aod_tevjets_v3_gidgt1_rt1cc1_p3_kWDW_v21_diag.root"

#
#hl_r3_24C = [
#    ["SRO_Data_Hist_sigma","",egres_r3_24_EG01_v12,"SRU"],
#    #["DRO_Data_Hist_sigma","",egres_r3_24_EG01_v12,"DRU"],
#    ["ZEE_Data_Hist_sigma","",egres_r3_24_EG01_v12,"ZEE"],
#]
#
#hl_r3_24d_part = [
#    ["SRO_Data_Hist_sigma","",egres_r3_23d_part_EG1_v12,"SRU"],
#    #["DRO_Data_Hist_sigma","",egres_r3_23d_part_EG1_v12,"DRU"],
#    ["ZEE_Data_Hist_sigma","",egres_r3_23d_part_EG1_v12,"ZEE"],
#]
#
#hl_r3_24d_part_cc = [
#    ["SRO_CC_Data_Hist_sigma","",egres_r3_23d_part_EG1_v12,"SRU"],
#    #["DRO_CC_Data_Hist_sigma","",egres_r3_23d_part_EG1_v12,"DRU"],
#    ["ZEE_CC_Data_Hist_sigma","",egres_r3_23d_part_EG1_v12,"ZEE"],
#]
#

hl = [
 
    #[ "gid1energy","",histfile3,"kOOT : True Ratio & True CC(3ns)"],
    #[ "gid1energy","",histfile12,"kOOT : True Ratio & True CC(12ns) "],
    #[ "gid23energy","",histfile3,"kOOT : False Ratio & False CC(3ns)"],
    #[ "gid23energy","",histfile12,"kOOT : False Ratio & False CC(12ns)"],
    #[ "gid2energy","",histfile3,"kOOT : False Ratio & True CC(3ns)"],
    #[ "gid2energy","",histfile12,"kOOT : False Ratio & True CC(12ns)"],
    [ "gidrfenergy","",histfile3,"kOOT : False Ratio"],
    [ "gidcfenergy","",histfile12,"kOOT : False CC(12ns)"],
    [ "gidcfenergy","",histfile3,"kOOT : False CC(3ns)"],
    [ "gidrtenergy","",histfile3,"kOOT : True Ratio"],
    [ "gidctenergy","",histfile12,"kOOT : True CC(12ns)"],
    [ "gidctenergy","",histfile3,"kOOT : True CC(3ns)"],
    #[ "gid3energy","",histfile3,"kOOT : True Ratio & False CC(3ns)"],
    #[ "gid3energy","",histfile12,"kOOT : True Ratio & False CC(12ns)"],
]

ptitle=[' 2024G TeVJet','','#splitline{EBEB}{gainID > 1}'] #{GT 106X_dataRun2_v28}'
y = [ 500000, 1 ]
#y = [ 1000, 1 ]
x = [ 0.0, 2000.0 ]
#l = [ 0.45,0.65,0.9,0.9 ]
#l = [ 0.45,0.65,0.9,0.9 ]
l = [ 0.65,0.65,0.9,0.9 ]
#t = [0.175,0.85,0.0,0.175,0.24]
t = [0.175,0.825,0.0,0.175,0.24]
outname = 'tr_hl_r3_24g_jetht_koot_false'
dostack(hl, outname, date, Ic_layout, ptitle,  y, x, l, t)


##ptitle=[' 2024C 379415-380238','','#splitline{EBEB}{5 GeV Cali}'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024C 379367-379543','','#splitline{EBEB}{ }'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024C 379434-380065','','#splitline{EBEB}{ }'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024C 380066-380238','','#splitline{EBEB}{ }'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024D 380306-380947','','#splitline{EBEB}{ }'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024D 380306-380947','','#splitline{EBEB}{5 GeV Cali}'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024CD 379415-380947','','#splitline{EBEB}{ }'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024CD 379415-380947','','#splitline{EBEB}{5 GeV Cali}'] #{GT 106X_dataRun2_v28}'
##ptitle=[' 2024E 381384','','#splitline{EBEB}{}'] #{GT 106X_dataRun2_v28}'
#ptitle=[' 2024D 370293-370580','','#splitline{Ratio timing}{EBEB CCGT}'] #{GT 106X_dataRun2_v28}'
#ptitle_cc=[' 2024D 370293-370580','','#splitline{CC timing}{EBEB CCGT}'] #{GT 106X_dataRun2_v28}'
#
#y = [ 0.4, 0.04 ]
#x = [ 50.0, 1200.0 ]
#l = [ 0.8,0.75,0.925,0.9 ]
#t = [0.175,0.44,0.0,0.175,0.24]
##outname = 'downloads/tr_hl_r3_24c_gold_v7'
##outname = 'downloads/tr_hl_r3_24e_test_v7'
#outname = 'downloads/tr_hl_r3_24d_part_ccgt_v7'
#outname_cc = 'downloads/tr_hl_r3_24d_part_ccgt_v7_cc'
#dostack(hl_r3_24d_part, outname, date, Ic_layout, ptitle,  y, x, l, t)
#dostack(hl_r3_24d_part_cc, outname_cc, date, Ic_layout, ptitle_cc,  y, x, l, t)

#
##    legend = TLegend(0.25,0.20,0.52,0.525); # bottom left
##    legend = TLegend(0.4,0.205,0.6,0.525);   # bottom middle
##    legend = TLegend(0.4,0.60,0.6,0.90);   # top middle
##    legend = TLegend(0.645,0.50,0.825,0.9);   # top mid right
##    legend = TLegend(0.605,0.50,0.945,0.9);   # top right very wide
##    legend = TLegend(0.705,0.50,0.945,0.9);   # top right wide 
##    legend = TLegend(0.745,0.50,0.925,0.9);   # top right
##    legend = TLegend(0.745,0.40,0.925,0.9);   # top right tall
##    legend = TLegend(0.650,0.375,0.925,0.875);   # top mid right wide
##    legend = TLegend(0.62,0.60,0.8,0.9);   # top right
##    legend = TLegend(0.65,0.60,0.9,0.90);   # top right large
#
##      g_2->GetYaxis()->SetMoreLogLabels();
##      g_2->GetYaxis()->SetNoExponent();
#

