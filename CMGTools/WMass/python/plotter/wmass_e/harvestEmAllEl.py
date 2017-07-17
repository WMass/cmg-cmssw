import CombineHarvester.CombineTools.ch as ch
import CombineHarvester.CombinePdfs.morphing as morphing
import ROOT
import glob, datetime, os

# import some parameters from wmass_parameters.py, they are also used by other scripts
from wmass_parameters import *

def harvestEm(subdir, mwrange, charge='both'):
    cmb = ch.CombineHarvester()
    
    # Read all the cards.
    # CH stores metadata about each object (Observation, Process, Systematic),
    # this is extracted from the card names with some regex
    for card in glob.glob(subdir+'/wenu_mass*.txt'): #.format(ch='' if charge == 'both' else 'plus*' if charge == 'plus' else 'minus*')):
        cmb.QuickParseDatacard(card, """wenu_mass(?<MASS>\d+)_$CHANNEL.card.txt""")
    
    # Need a unqiue bin name for each plus/minus,pt and eta combination
    # We extracted this part of the datacard name into the channel variable above,
    # so can just copy it and override the specific bin name that was in all the cards
    cmb.ForEachObj(lambda obj: obj.set_bin(obj.channel()))
 
    # We'll have three copies of the observation, one for each mass point.
    # Filter all but one copy.
    #cmb.FilterObs(lambda obj: obj.mass() != '15')
    cmb.FilterObs(lambda obj: obj.mass() != str(mass_id_central))
    
    # Create workspace to hold the morphing pdfs and the mass
    w = ROOT.RooWorkspace('morph', 'morph')
    mass = w.factory('mw[{mwrange}]'.format(mwrange=mwrange))
    
    # BuildRooMorphing will dump a load of debug plots here
    debug = ROOT.TFile(subdir+'/debug.root', 'RECREATE')
    
    # Run for each bin,process combination (only for signal!)
    for b in cmb.bin_set():
        for p in cmb.cp().bin([b]).signals().process_set():
            morphing.BuildRooMorphing(w, cmb, b, p, mass, verbose=True, file=debug)
    
    # Just to be safe
    mass.setConstant(True)
    
    # Now the workspace is copied into the CH instance and the pdfs attached to the processes
    # (this relies on us knowing that BuildRooMorphing will name the pdfs in a particular way)
    cmb.AddWorkspace(w, True)
    cmb.cp().process(['W']).ExtractPdfs(cmb, 'morph', '$BIN_$PROCESS_morph', '')
    
    # Adjust the rateParams a bit - we currently have three for each bin (one for each mass),
    # but we only want one. Easiest to drop the existing ones completely and create new ones
    cmb.syst_type(['rateParam'], False)
    cmb.cp().process(['W']).AddSyst(cmb, 'norm_$BIN', 'rateParam', ch.SystMap()(1.00))
    
    # Have to set the range by hand
    for sys in cmb.cp().syst_type(['rateParam']).syst_name_set():
        cmb.GetParameter(sys).set_range(0.5, 1.5)
    
    # Print the contents of the model
    cmb.PrintAll()
    
    # Write out the cards, one per bin
    outdir=subdir+'/wenu_cards_morphed_{charge}'.format(charge=charge)
    writer = ch.CardWriter('$TAG/$BIN.txt', '$TAG/shapes.root')
    writer.SetVerbosity(1)
    writer.WriteCards(outdir, cmb)

date = datetime.date.today().isoformat()
date+='_charges'

card_dir = 'cards/lepPtOptim_12massVar_ptLow30/'
subdirs = [x[0] for x in os.walk(card_dir)]

#mwrange='0,30'
# values from wmass_parameters.py
mwrange='%d,%d' % (mass_id_down,mass_id_up)
npoints = n_mass_id
central = mass_id_central

runHarvest = False
runBatch   = False
justHadd   = False
combineCards = False
runFit = True

input_dcs_alleta = ""
workspaces = []
for isub, subdir in enumerate(subdirs):
    if subdir == subdirs[0]: continue
    if 'wenu_cards_morphed' in subdir: continue
    name = subdir.split('/')[-1]
    if not 'eta_' in name: continue

    print '--------------------------------------------------------------------'
    print '- running for {mode} -----------------------------------------------'.format(mode=name)
    print '- in subdirectory {subdir} -----------------------------------------'.format(subdir=subdir)
    print '--------------------------------------------------------------------'
    #if name == 'full_3d': continue

    if runHarvest: 
        ## run the combine harvester which combines all the datacards etc.
        harvest(subdir,mwrange)

    target_dc = '{subdir}/wenu_cards_morphed_both/morphed_datacard_channel.txt'.format(subdir=subdir)
    target_ws = target_dc.replace('txt','root')
    workspaces.append(target_ws)

    if combineCards:
        ## running combineCards to make the combined plus+minus datacard
        if os.path.isfile(target_dc):
            print 'removing existing combined datacard first!'
            os.system('rm {dc}'.format(dc=target_dc) )
        dcs = os.listdir(subdir+"/wenu_cards_morphed_both/")
        input_dcs=" ".join(["%s=%s" % (os.path.splitext(dc)[0],subdir+"/wenu_cards_morphed_both/"+dc) for dc in dcs if "txt" in dc])
        input_dcs_alleta += " "+input_dcs

        print 'running combineCards.py'
        combineCardsCmd = 'combineCards.py {dcs} >& {target_dc}'.format(dcs=input_dcs, target_dc=target_dc)
        print combineCardsCmd
        ## run combineCards and make the workspace
        os.system(combineCardsCmd )
        print 'running text2workspace'
        t2wCmd = 'text2workspace.py {target_dc} '.format(subdir=subdir, target_dc=target_dc)
        print t2wCmd
        os.system(t2wCmd)

comb_dir = card_dir+'/comb'
if not os.path.exists(comb_dir): os.mkdir(comb_dir)
comb_dc = comb_dir+"/morphed_datacard_comb.txt"
comb_ws = comb_dc.replace('txt','root')
workspaces.append(comb_ws)
        
if combineCards:
    if os.path.isfile(comb_dc):
        print 'removing existing combined datacard first!'
        os.system('rm {dc}'.format(dc=comb_dc) )

    print 'running combineCards.py'
    combineCardsCmd = 'combineCards.py {dcs} >& {target_dc}'.format(dcs=input_dcs_alleta, target_dc=comb_dc)
    print combineCardsCmd
    ## run combineCards and make the workspace
    os.system(combineCardsCmd)
    print 'running text2workspace'
    os.system('text2workspace.py %s' % comb_dc)

if runFit:
    for ws in workspaces:
        print "===> RUN FIT FOR WORKSPACE: ",ws
        ## constructing the command
        combine_base  = 'combine -t -1 -M MultiDimFit --setPhysicsModelParameters mw={central},r=1 --setPhysicsModelParameterRanges mw={mwrange} '.format(central=central,mwrange=mwrange)
        combine_base += ' --redefineSignalPOIs=mw --algo grid --points {npoints} {target_ws} '.format(npoints=npoints, target_ws=ws)
        
        saveNuisances = ''
        saveNuisances += ' --saveSpecifiedNuis {vs} '.format(vs=','.join('CMS_We_pdf'+str(i) for i in range(1,27)))
        
        run_combine_pdfUnc = combine_base + ' -n {date}_{name} {sn} '.format(date=date,name=name,sn=saveNuisances) 
        run_combine_noPdf  = combine_base + ' -n {date}_{name}_noPDFUncertainty --freezeNuisanceGroups pdfUncertainties '.format(date=date,name=name)
        
        if runBatch:
            run_combine_pdfUnc += ' --job-mode lxbatch --split-points 10 --sub-opts="-q 8nh" --task-name {name}                  '.format(name=name)
            run_combine_noPdf  += ' --job-mode lxbatch --split-points 10 --sub-opts="-q 8nh" --task-name {name}_noPDFUncertainty '.format(name=name)
            run_combine_pdfUnc  = 'combineTool.py ' + ' '.join(run_combine_pdfUnc.split()[1:])
            run_combine_noPdf   = 'combineTool.py ' + ' '.join(run_combine_noPdf .split()[1:])
        
        
        ## running combine once with the systematics and once without
        print '-- running combine command ------------------------------'
        print '---     with uncertainties: -----------------------------'
        print run_combine_pdfUnc
        
        os.system(run_combine_pdfUnc)
        
        print '---     without uncertainties: --------------------------'
        print run_combine_noPdf
        os.system(run_combine_noPdf )
        
        impactBase = 'combineTool.py -M Impacts -n {date}_{name} -d {target_ws} -m {mass} '.format(mass=name[-1],date=date,name=name, target_ws=ws)
        impactBase += ' --setPhysicsModelParameters mw={central},r=1  --redefineSignalPOIs=mw --setPhysicsModelParameterRanges mw={mwrange} -t -1 '.format(central=central,mwrange=mwrange)
        impactInitial = impactBase+'  --robustFit 1 --doInitialFit '
        impactFits    = impactBase+'  --robustFit 1 --doFits '
        impactJSON    = impactBase+'  -o impacts_{name}.json '.format(name=name)
        impactPlot    = 'plotImpacts.py -i impacts_{name}.json -o impacts_{name} --transparent'.format(name=name)
        
        os.system(impactInitial)
        os.system(impactFits   )
        os.system(impactJSON   )
        os.system(impactPlot   )

## combineTool.py -M Impacts -d morphed_datacard_plusminus.root -m 999 --setPhysicsModelParameters mw=19,r=1  --redefineSignalPOIs=mw  --setPhysicsModelParameterRanges mw=0,30 -t -1 --robustFit 1 --doInitialFit
## combineTool.py -M Impacts -d morphed_datacard_plusminus.root -m 999 --setPhysicsModelParameters mw=19,r=1  --redefineSignalPOIs=mw  --setPhysicsModelParameterRanges mw=0,30 -t -1 --robustFit 1 --doFits
## combineTool.py -M Impacts -d morphed_datacard_plusminus.root -m 999 --setPhysicsModelParameters mw=19,r=1  --redefineSignalPOIs=mw  --setPhysicsModelParameterRanges mw=0,30 -t -1 -o impacts.json
## plotImpacts.py -i impacts.json -o impacts


