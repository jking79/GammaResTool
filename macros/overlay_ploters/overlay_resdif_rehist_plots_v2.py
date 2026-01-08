#  jack w king 3

#from numpy import array
from ROOT import *
from math import *
#from helper import *  
from os import system as call
import time

from jwk_tdr_style_py import *

def dostack( hist_list, outname, date, layout, ptitle, y, x, l, t ):

    dofit = True
    #dofit = False
    paramn = []
    parnerror = []
    paramc = []
    parcerror = []
    params = []
    parserror = []

    first = True
    f1 = []
    f2 = []
    f1b = []
    f2b = []
    h1 = []
    hfit = []

    setTDRStyle()
    gStyle.SetPadTickY(1)
    gStyle.SetPadTickX(1)

    #gStyle.SetPaperSize(20,26);
    gStyle.SetPadTopMargin(0.09);
    gStyle.SetPadRightMargin(0.05);
    gStyle.SetPadBottomMargin(0.15);
    gStyle.SetPadLeftMargin(0.15);

    if layout['logx'] : gStyle.SetOptLogx(1)
    if layout['logy'] : gStyle.SetOptLogy(1)
    c1 = TCanvas( 'c1', 'canvas' , 200, 10, 700, 500 )
    c1.cd()
    #if layout['logx'] : c1.SetLogx()
    #if layout['logy'] : c1.SetLogy()
    c1.SetGridx(1)
    c1.SetGridy(1)
    c1.Update()
    c1.Draw()

    legend = TLegend(l[0],l[1],l[2],l[3]);
#    legend.SetHeader(layout['legtitle'], '')
    legend.SetName('legend')
    gStyle.SetLegendFont(42)
    gStyle.SetLegendTextSize(0.03)
    legend.SetBorderSize(0)
    #legend.SetLineColor(kBlack)

    lat = TLatex() 
    lat.SetNDC()

    hMax = y[0]
    hMin = y[1]

    for n in range(0,3):

        #ns=str(n)
        hfit.append(TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] ) )',75,750,2))
        #hfit = TF1('hfits','sqrt( ( ([0]*[0])/(x*x) )+( 2*[1]*[1] )+( ([2]*[2])/(x) ) )',50,450,3)
        hfit[n].SetParName(0,'N')
        hfit[n].SetParameter(0,40.0)
        #hfit.SetParLimits(0,0,50)
        hfit[n].SetParLimits(0,0.0,100.0)
        #hfit.SetParameter(0,5)
        #hfit.SetParLimits(0,0,10)
        hfit[n].SetParName(1,'C')
        hfit[n].SetParameter(1,0.1)
        hfit[n].SetParLimits(1,0.0,1.0)
        #hfit.SetParLimits(1,0.01,10.0)
        #hfit.SetParLimits(1,0.02,1.0)
        #hfit.SetParameter(1,0.05)
        #hfit.SetParLimits(1,0.001,1.0)
        #hfit.SetParName(2,'S')
        #hfit.SetParameter(2,5.0)
        #hfit.SetParLimits(2,0.0,25.0)

    mg = TMultiGraph();
    first = True;

    for histname, histname2, infileb1, infile1, lego in hist_list :
    
        f1.append(TFile.Open(infileb1))
        #if tree == '' : hist = histname
        #else : hist = tree+'/'+histname

        f2.append(TFile.Open(infile1))
        #if tree == '' : hist = histname2
        #else : hist = tree+'/'+histname2
 
        orighist1 = f1[0].Get(histname)
        orighist2 = f2[0].Get(histname2)
        #htitle = 'hist' + str(n)
        lenmybins = int(orighist1.GetNbinsX())
        print( "length: ", lenmybins )
        h1.append(TGraphErrors(orighist1))
        h1.append(TGraphErrors(orighist2))
        h1.append(TGraphErrors(orighist2))
        #h1[0] = TGraphErrors(lenmybins)
        #h1[1] = TGraphErrors(lenmybins)
        #h1[2] = TGraphErrors(lenmybins)
        #norm1 = orighist1.Integral()
        #norm2 = orighist2.Integral()
        #if norm1 == 0 : norm1 = 1
        #if norm2 == 0 : norm2 = 1
        #print( "Norms : ", norm1, norm2 )

        for bn in range( 0, lenmybins):

            binmid1 = float(orighist1.GetBinCenter(bn))
            #binmid2 = float(orighist2.GetBinCenter(bn))
            obinval1 = float(orighist1.GetBinContent(bn))
            obinval2 = float(orighist2.GetBinContent(bn))
            binerr1 = float(orighist1.GetBinError(bn))
            binerr2 = float(orighist2.GetBinError(bn))
            binwidth1 = float(orighist1.GetBinWidth(bn))
            #binwidth2 = float(orighist2.GetBinWidth(bn))
            diff = obinval1*obinval1 - obinval2*obinval2
            if diff < 0 : diff = -1*diff
            valdiff = sqrt( diff )
            print( binmid1, obinval1, obinval2, binerr1, binerr2, valdiff )
            errdif = valdiff
            if errdif == 0 : errdif = 1
            err = sqrt( binerr1*binerr1*obinval1*obinval1 + binerr2*binerr2*obinval2*obinval2 )/errdif
            #h1[0].SetPoint(bn,binmid1,obinval1)
            #h1[1].SetPoint(bn,binmid1,obinval2)
            h1[2].SetPoint(bn,binmid1,valdiff)
            #h1[n+1].SetPoint(bn,binmid2,obinval2)
            #if obinval1 == 0 : obinval1 = 1;
            #if obinval2 == 0 : obinval2 = 1; 
            widtherr1 = 0 #binwidth1/(2*sqrt(obinval1))
            #widtherr2 = 0 #binwidth2/(2*sqrt(obinval2))
            #h1[0].SetPointError(bn,widtherr1,binerr1)
            #h1[1].SetPointError(bn,widtherr1,binerr2)
            h1[2].SetPointError(bn,widtherr1,err)
            #h1[n+1].SetPointError(bn,widtherr2,binerr2)
            #h1[n].SetPoint(bn,sigpass,bkgpass)
            #print( " set point : ", bn, sigpass, bkgpass,  )

#---------------------------------------------------------------------

    for n in range(0,3):

        h1[n].UseCurrentStyle()
        h1[n].SetMarkerStyle(n+25)
        #k = [kMagenta+2,kGreen+2,kYellow+1,kBlue+2,kRed+2,kAzure+4,kBlack,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2]
        k = [kMagenta+2,kGreen+2,kBlue+2,kRed+2,kAzure+4,kViolet+7,kOrange+7,kGreen+3,kRed+4,kBlue+4,kGreen+2,kAzure+4,kMagenta+2,kGreen+2,kBlack]
        #k = [kBlue+4,kBlue+1,kGreen+4,kYellow+3,kAzure+4,kViolet+7,kOrange+7,kGreen+3]
        #k = [kSpring-7,kSpring+3,kAzure+3,kAzure-7]
        #k = [kBlack]
        #k = [kGray+1,kGray+2,kGray+3,kBlack]
        h1[n].SetLineColor(k[n])
        h1[n].SetMarkerColor(k[n])
        if dofit :
            hfit[n].SetLineColor(k[n])
            hfit[n].SetLineWidth(2)
        msz = 0.8
        if( n == 1 ) : h1[n].SetMarkerSize(msz+0.3)
        elif( n == 2 ) : h1[n].SetMarkerSize(msz+0.5)
        else : h1[n].SetMarkerSize(msz)
        h1[n].SetMarkerSize(msz)

        if( first ) :
                if dofit : h1[n].Fit(hfit[n],'RE')
                first = False
        else :
                if dofit : h1[n].Fit(hfit[n],'RE+')

        if dofit :
                 paramn.append(str(abs(hfit[n].GetParameter(0))))
                 paramc.append(str(abs(hfit[n].GetParameter(1))))
                 #params.append(str(abs(hfit.GetParameter(2))))
                 pne = hfit[n].GetParError(0)
                 pce = hfit[n].GetParError(1)
                 #pse = hfit.GetParError(2)
                 print('Fit info',paramn[n],pne,paramc[n],pce)
                 #print('Fit info',paramn[n],pne,paramc[n],pce,params[n],pse)
                 if pne < 0.01 : pne = 0.01
                 if pce < 0.0001 : pce = 0.0001
                 parnerror.append(str(pne))
                 parcerror.append(str(pce))
                 #parserror.append(str(pse))

        #if n == 0 : lego = "Online"
        #if n == 1 : lego = "+Offline"
        #if n == 2 : lego = "Difference"
        if n == 0 : lego = "Data (DEG 17E)"
        if n == 1 : lego = "MC (DY)"
        if n == 2 : lego = "Diffrence"
        legend.AddEntry(h1[n],lego,'epl');
        n += 1

        #End of loop
    
    for n in range(0,3):
        mg.Add(h1[n])

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
    mg.SetMinimum(y[0])
    mg.SetMaximum(y[1])
#   mg.GetXaxis().SetRangeUser(200.0,1100.0)
    mg.GetXaxis().SetRangeUser(x[0],x[1])
#    if layout['logx'] : mg.GetXaxis().SetMoreLogLabels()
#    if layout['logy'] : mg.GetYaxis().SetMoreLogLabels()

    if lego != 'none' : legend.Draw('same')  #   legend inclusion switch
    gPad.Modified()

    #lat_cms = '#bf{CMS} #it{Preliminary}' + ptitle[0]
    lat_cms = '#bf{CMS} #it{Work in Progress}' + ptitle[0]
    #lat_title = ptitle[1]+' (13 TeV)'
    lat_title = ptitle[1]
    lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (#frac{N}{A_{eff}/#sigma_{n}})^{2} + #frac{S^{2}}{A_{eff}/#sigma_{n}} + 2C^{2}'
    #lat_form = '#sigma^{2}_{i} = (N/Eeff)^{2} + 2C^{2}'
    lat.SetTextSize(0.05);
    lat.SetTextFont(42);
    lat.DrawLatex(0.15,0.9325,lat_cms);
    lat.DrawLatex((0.82-t[2]),0.9325,lat_title);
    lat.SetTextSize(0.04);
    lat.DrawLatex(t[0],t[1],ptitle[2]);
   
    if dofit :
    #if False :
        lat.SetTextSize(0.03);
        lat.DrawLatex(t[3],t[4]+.08,lat_form);
        for l in range(0,3):
            lat_param = '#color['+str(k[l])+']{'
            lat_param = lat_param + 'N : '+paramn[l][0:4]+' #pm '+parnerror[l][0:3]+' [ns]   '
            #lat_param = lat_param + 'S : '+params[l][0:3]+' #pm '+parserror[l][0:3]+' [ns]   '
            lat_param = lat_param + 'C : '+paramc[l][0:6]+' #pm '+parcerror[l][0:6]+' [ns]}'
            lat.SetTextSize(0.03);
            lat.DrawLatex(t[3],t[4]-l*.035,lat_param);
 
    if layout['logx'] : c1.SetLogx()
    c1.Modified()
    c1.Update()

    #c1.Print( outname + '_' + date + '.pdf' )
    c1.Print( outname + '_' + date + '.png' )
    #c1.SaveAs( outname + '_' + date + '.root' )
    #c1.SaveAs( outname + '_' + date + '.C' )
    #c1.Show()
    c1.Close()

#from overlay_hist_defs_v3 import *

date = time.strftime('%d%m%Y') +  time.strftime('%H%M%S')

legtitle = ''
#legtitle = 'KuStc'
#legtitle = 'KuNotStc'

rtitle = 'Run2 AODSIM'
#Ic_legtitle = ''
xtitle = 'A_{eff}/#sigma_{n}'
#xtitle = '#Delta_{Run}'
#xtitle = 'GeV'
#xtitle = '[ns]'
#ytitle = 'Ave Xtal Time [ns]'
#ytitle = '#sigma(Adjusted pCalo time) [ns]'
#ytitle = '#mu(t_{1}-t_{2}) [ns]'
#ytitle = '#occupancy(t_{1}-t_{2}) '
#xtitle = 'HcalTowerSumEtBcConeDR04 [GeV]'
#ytitle = 'a.u.'
ytitle = '#sigma_{t_{1}-t_{2}} [ns]'
#ytitle = '#Delta(#sigma_{DEG},#sigma_{DY}) [ns]'
htitle = ''
#islogx = True
islogx = False
#islogy = True
islogy = False

#---------------------------------------------------------------
#hl_mc_fms_loc = [
#     #['hist_name","tree_name",hist_file_location","legend_name"],
#     #["Data_sigma","",mc_full_loc+pcal+lstfr,"Full"],
#     #["Data_sigma","",mc_multi_loc+pcal+lstfr,"Multi"],
#     #["Data_sigma","",mc_single_loc+pcal+lstfr,"Single"],
#]

l = [ 0.7,0.7,0.925,0.9 ] # legend position top right
#l = [ 0.2,0.55,0.425,0.75 ] # legend position top left
#l = [ 0.25,0.20,0.52,0.525 ] # legend position bottom left
t = [0.3,0.84,0.0,0.175,0.25] # titles position

rhname1 = 'ZEE_Data_Hist_sigma'

#degdir = 'egres_DEGPD_AOD_Run2017F_Cali_FilteredFit_Full_14011_v10_plots/'
#degres = 'egres_DEGPD_AOD_Run2017F_v2_Cali_FilteredFit_Full_14011_v12_resplots.root'
#dydir = 'egres_DY1JetsToLL_AODSIM_R17_Cali_FilterFit_14011_v10_plots/'
#dyres = 'egres_DY1JetsToLL_AODSIM_R17_FilterFit_14011_v12_resplots.root'

#degdir = 'egres_DEGPD_AOD_Run2017F_Cali_FilteredFit_ZBins_Full_14011_v10_plots/'
#degres = 'egres_DEGPD_AOD_Run2017F_v2_Cali_FilteredFit_ZBins_Full_14011_v12_resplots.root'
degdir = 'egres_DEGPD_AOD_Run2017EF_Cali_FilteredFit_ZBins_Full_14011_v11_plots/'
degres = 'egres_DEGPD_AOD_Run2017EF_v2_Cali_FilteredFit_ZBins_Full_14011_v12_resplots.root'
dydir = 'egres_DY1JetsToLL_AODSIM_R17_Cali_FilterFit_ZBins_14011_v10_plots/'
dyres = 'egres_DY1JetsToLL_AODSIM_R17_FilterFit_ZBins_14011_v12_resplots.root'

egrzcali = 'ecal_config/EG_EOY_MINI_304476_X_ZEE_Data_Hist_resfit.root'
egrznocali = 'ecal_config/EG_EOY_MINI_304476_X_ZEE_Data_Hist_NoCali_resfit.root'

rhname1 = 'EG_EOY_MINI_304476_X_ZEE_Data_Hist_sigma'
rhname2 = 'EG_EOY_MINI_304476_X_ZEE_Data_Hist_NoCali_sigma'

resdir = 'res_files/'
degdir = 'egres_DEGPD_AOD_Run2017EF_Cali_FilteredFit_ZBins_Full_14011_v11_plots/'
degres = 'egres_DEGPD_AOD_Run2017EF_v2_Cali_FilteredFit_ZBins_Full_14011_v12_resplots.root'
dydir = 'egres_DY1JetsToLL_AODSIM_R17_Cali_FilterFit_ZBins_14011_v10_plots/'
dyres = 'egres_DY1JetsToLL_AODSIM_R17_FilterFit_ZBins_14011_v12_resplots.root'

loc1 = resdir + degdir + degres
loc2 = resdir + dydir + dyres

rhname3 = 'ZEE_Data_Hist_sigma'

x = [ 60, 750 ]
y = [ 0, 0.5 ]

title = ''
#title = '#splitline{EBEB Z #rightarrow ee}{}'
ptcut = ' '
pt = ''
#isocut = '#splitline{#splitline{hadTowOverEM < 0.02}{ecalRHSumEtConeDR04 < 9.0}}{trkSumPtSolidConeDR04 < 6.0}'
isocut = '#splitline{EBEB}{Z#rightarrow e^{+}e^{-}}'
#ptlist = [ '20','30','100']
#outfileend = '_pt'

#for rhname1, rhname2, rhname3, rhname4, pt in zip( rhname1l, rhname2l, rhname3l, rhname4l, ptlist ) :

if True :

  ptitle=[ title, ptcut+pt, isocut ]
  #outname = 'llpa_egmres_calidiff_zee'
  outname = 'llpa_egmres_smeardiff_zee'
  layout = { 'xtitle' : xtitle, 'ytitle' : ytitle, 'title' : htitle, 'logx' : islogx, 'logy' : islogy, 'legtitle' : legtitle }
  
  inhistlist = [
  
  #    [ rhname1, "", base1-sig, base2-bkgrd, rfgmsbroc1-sig, rfgmsbroc2-bkgrd, rhname1, cut1, usebase1 ],
      #[ rhname1, rhname2, egrzcali, egrznocali, '' ],
      [ rhname3, rhname3, loc1, loc2, '' ],  

  ]
  
  dostack(inhistlist, outname, date, layout, ptitle,  y, x, l, t)
  

#ptitle=[' 2022 IOV5 359421-360089','','#splitline{EBEB}{CC Ave RH Time by Channel}'] #{GT 106X_dataRun2_v28}'
#y = [ 4.5, 0.5 ]
#x = [ 200.0, 700.0 ]
#l = [ 0.7,0.65,0.925,0.9 ]
#t = [0.2,0.825,0.0,0.175,0.225]
#outname = 'downloads/tr_hl_r3_iov5cali_v7'
#dostack(hl_r3_iov5cali_v7, outname, date, Ic_layout, ptitle,  y, x, l, t)

#    legend = TLegend(0.25,0.20,0.52,0.525); # bottom left
#    legend = TLegend(0.4,0.205,0.6,0.525);   # bottom middle
#    legend = TLegend(0.4,0.60,0.6,0.90);   # top middle
#    legend = TLegend(0.645,0.50,0.825,0.9);   # top mid right
#    legend = TLegend(0.605,0.50,0.945,0.9);   # top right very wide
#    legend = TLegend(0.705,0.50,0.945,0.9);   # top right wide 
#    legend = TLegend(0.745,0.50,0.925,0.9);   # top right
#    legend = TLegend(0.745,0.40,0.925,0.9);   # top right tall
#    legend = TLegend(0.650,0.375,0.925,0.875);   # top mid right wide
#    legend = TLegend(0.62,0.60,0.8,0.9);   # top right
#    legend = TLegend(0.65,0.60,0.9,0.90);   # top right large

#      g_2->GetYaxis()->SetMoreLogLabels();
#      g_2->GetYaxis()->SetNoExponent();


