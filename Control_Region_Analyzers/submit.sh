# Compile
./rootcom postAnalyzer_ZnnG_mc_znng ZnnG_mc_znng
./rootcom postAnalyzer_ZnnG_mc_wg ZnnG_mc_wg
./rootcom postAnalyzer_ZnnG_mc ZnnG_mc
./rootcom postAnalyzer_ZnnG_data ZnnG_data
./rootcom postAnalyzer_ZnnG_wenu ZnnG_wenu
./rootcom postAnalyzer_ZnnG_qcd ZnnG_qcd
./rootcom postAnalyzer_ZnnG_bhalo ZnnG_bhalo

./rootcom postAnalyzer_ZllG_mc_znng ZllG_mc_znng
./rootcom postAnalyzer_ZllG_mc_wg ZllG_mc_wg
./rootcom postAnalyzer_ZllG_mc ZllG_mc
./rootcom postAnalyzer_ZllG_data ZllG_data
./rootcom postAnalyzer_ZllG_wenu ZllG_wenu

./rootcom postAnalyzer_WenG_mc_znng WenG_mc_znng
./rootcom postAnalyzer_WenG_mc_wg WenG_mc_wg
./rootcom postAnalyzer_WenG_mc WenG_mc
./rootcom postAnalyzer_WenG_data WenG_data
./rootcom postAnalyzer_WenG_wenu WenG_wenu
./rootcom postAnalyzer_WenG_qcd WenG_qcd

./rootcom postAnalyzer_WmnG_mc_znng WmnG_mc_znng
./rootcom postAnalyzer_WmnG_mc_wg WmnG_mc_wg
./rootcom postAnalyzer_WmnG_mc WmnG_mc
./rootcom postAnalyzer_WmnG_data WmnG_data
./rootcom postAnalyzer_WmnG_wenu WmnG_wenu
./rootcom postAnalyzer_WmnG_qcd WmnG_qcd


# Run
# ewk_corr.root does not play nicely with condor, so run anything that needs it locally
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  ZnnG_data0000b.root -1 1000 ZnnG_data0000b
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  ZnnG_data0001b.root -1 1000 ZnnG_data0001b
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  ZnnG_data0002b.root -1 1000 ZnnG_data0002b
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  ZnnG_data0003b.root -1 1000 ZnnG_data0003b
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  ZnnG_data0000c.root -1 1000 ZnnG_data0000c
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  ZnnG_data0001c.root -1 1000 ZnnG_data0001c
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  ZnnG_data0000d.root -1 1000 ZnnG_data0000d
./MakeCondorFiles.csh ZnnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  ZnnG_data0001d.root -1 1000 ZnnG_data0001d
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  ZnnG_wenu0000b.root -1 1000 ZnnG_wenu0000b
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  ZnnG_wenu0001b.root -1 1000 ZnnG_wenu0001b
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  ZnnG_wenu0002b.root -1 1000 ZnnG_wenu0002b
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  ZnnG_wenu0003b.root -1 1000 ZnnG_wenu0003b
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  ZnnG_wenu0000c.root -1 1000 ZnnG_wenu0000c
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  ZnnG_wenu0001c.root -1 1000 ZnnG_wenu0001c
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  ZnnG_wenu0000d.root -1 1000 ZnnG_wenu0000d
./MakeCondorFiles.csh ZnnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  ZnnG_wenu0001d.root -1 1000 ZnnG_wenu0001d
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  ZnnG_qcd0000b.root -1 1000 ZnnG_qcd0000b
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  ZnnG_qcd0001b.root -1 1000 ZnnG_qcd0001b
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  ZnnG_qcd0002b.root -1 1000 ZnnG_qcd0002b
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  ZnnG_qcd0003b.root -1 1000 ZnnG_qcd0003b
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  ZnnG_qcd0000c.root -1 1000 ZnnG_qcd0000c
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  ZnnG_qcd0001c.root -1 1000 ZnnG_qcd0001c
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  ZnnG_qcd0000d.root -1 1000 ZnnG_qcd0000d
./MakeCondorFiles.csh ZnnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  ZnnG_qcd0001d.root -1 1000 ZnnG_qcd0001d

nohup ./ZnnG_mc_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/ ZnnG_ZNuNuGJets.root -1 1000 >& ZnnG_ZNuNuGJets.txt &
nohup ./ZnnG_mc_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ ZnnG_WGJets.root -1 1000 >& ZnnG_WGJets.txt &
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/160621_133030/0000/  ZnnG_GJets_HT-40To100.root -1 1000 ZnnG_GJets_HT-40To100
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/160621_133042/0000/  ZnnG_GJets_HT-100To200.root -1 1000 ZnnG_GJets_HT-100To200
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/160621_133055/0000/  ZnnG_GJets_HT-200To400.root -1 1000 ZnnG_GJets_HT-200To400
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/160621_133110/0000/  ZnnG_GJets_HT-400To600.root -1 1000 ZnnG_GJets_HT-400To600
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/160621_133123/0000/  ZnnG_GJets_HT-600ToInf.root -1 1000 ZnnG_GJets_HT-600ToInf
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WToMuNu_M-100_TuneCUETP8M1_13TeV-pythia8/crab_WToMuNu/160621_133001/0000/  ZnnG_WToMuNu.root -1 1000 ZnnG_WToMuNu
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WToTauNu_M-100_TuneCUETP8M1_13TeV-pythia8-tauola/crab_WToTauNu/160621_133014/0000/  ZnnG_WToTauNu.root -1 1000 ZnnG_WToTauNu
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/  ZnnG_ZLLGJets.root -1 1000 ZnnG_ZLLGJets
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/160621_133135/0000/  ZnnG_TTGJets.root -1 1000 ZnnG_TTGJets
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets_MGG-80toInf_amcatnlo/160825_000327/0000/  ZnnG_Diphoton.root -1 1000 ZnnG_Diphoton
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TGJets/160823_154537/0000/  ZnnG_TGJets.root -1 1000 ZnnG_TGJets
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_WWG/160630_154410/0000/  ZnnG_WWG.root -1 1000 ZnnG_WWG
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/160630_154436/0000/  ZnnG_ZZ.root -1 1000 ZnnG_ZZ
./MakeCondorFiles.csh ZnnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/160630_154423/0000/  ZnnG_WZ.root -1 1000 ZnnG_WZ

./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  ZllG_data0000b.root -1 1000 ZllG_data0000b
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  ZllG_data0001b.root -1 1000 ZllG_data0001b
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  ZllG_data0002b.root -1 1000 ZllG_data0002b
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  ZllG_data0003b.root -1 1000 ZllG_data0003b
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  ZllG_data0000c.root -1 1000 ZllG_data0000c
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  ZllG_data0001c.root -1 1000 ZllG_data0001c
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  ZllG_data0000d.root -1 1000 ZllG_data0000d
./MakeCondorFiles.csh ZllG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  ZllG_data0001d.root -1 1000 ZllG_data0001d
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  ZllG_wenu0000b.root -1 1000 ZllG_wenu0000b
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  ZllG_wenu0001b.root -1 1000 ZllG_wenu0001b
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  ZllG_wenu0002b.root -1 1000 ZllG_wenu0002b
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  ZllG_wenu0003b.root -1 1000 ZllG_wenu0003b
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  ZllG_wenu0000c.root -1 1000 ZllG_wenu0000c
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  ZllG_wenu0001c.root -1 1000 ZllG_wenu0001c
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  ZllG_wenu0000d.root -1 1000 ZllG_wenu0000d
./MakeCondorFiles.csh ZllG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  ZllG_wenu0001d.root -1 1000 ZllG_wenu0001d

nohup ./ZllG_mc_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ ZllG_WGJets.root -1 1000 >& ZllG_WGJets.txt &
nohup ./ZllG_mc_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ ZllG_ZLLGJets.root -1 1000 >& ZllG_ZLLGJets.txt &
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/160621_133135/0000/  ZllG_TTGJets.root -1 1000 ZllG_TTGJets
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-100to200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-100to200/160621_133238/0000/  ZllG_DYJetsToLL_HT-100to200.root -1 1000 ZllG_DYJetsToLL_HT-100to200
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-200to400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-200to400/160621_133251/0000/  ZllG_DYJetsToLL_HT-200to400.root -1 1000 ZllG_DYJetsToLL_HT-200to400
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-400to600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-400to600/160621_133304/0000/  ZllG_DYJetsToLL_HT-400to600.root -1 1000 ZllG_DYJetsToLL_HT-400to600
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DYJetsToLL_M-50_HT-600toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_DYJetsToLL_M-50_HT-600toInf/160621_133316/0000/  ZllG_DYJetsToLL_HT-600toInf.root -1 1000 ZllG_DYJetsToLL_HT-600toInf
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/160630_154423/0000/  ZllG_WZ.root -1 1000 ZllG_WZ
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/160630_154436/0000/  ZllG_ZZ.root -1 1000 ZllG_ZZ
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_WWG/160630_154410/0000/  ZllG_WWG.root -1 1000 ZllG_WWG
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/160621_133030/0000/  ZllG_GJets_HT-40To100.root -1 1000 ZllG_GJets_HT-40To100
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/160621_133042/0000/  ZllG_GJets_HT-100To200.root -1 1000 ZllG_GJets_HT-100To200
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/160621_133055/0000/  ZllG_GJets_HT-200To400.root -1 1000 ZllG_GJets_HT-200To400
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/160621_133110/0000/  ZllG_GJets_HT-400To600.root -1 1000 ZllG_GJets_HT-400To600
./MakeCondorFiles.csh ZllG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/160621_133123/0000/  ZllG_GJets_HT-600ToInf.root -1 1000 ZllG_GJets_HT-600ToInf

./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  WenG_data0000b.root -1 1000 WenG_data0000b
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  WenG_data0001b.root -1 1000 WenG_data0001b
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  WenG_data0002b.root -1 1000 WenG_data0002b
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  WenG_data0003b.root -1 1000 WenG_data0003b
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  WenG_data0000c.root -1 1000 WenG_data0000c
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  WenG_data0001c.root -1 1000 WenG_data0001c
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  WenG_data0000d.root -1 1000 WenG_data0000d
./MakeCondorFiles.csh WenG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  WenG_data0001d.root -1 1000 WenG_data0001d
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  WenG_wenu0000b.root -1 1000 WenG_wenu0000b
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  WenG_wenu0001b.root -1 1000 WenG_wenu0001b
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  WenG_wenu0002b.root -1 1000 WenG_wenu0002b
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  WenG_wenu0003b.root -1 1000 WenG_wenu0003b
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  WenG_wenu0000c.root -1 1000 WenG_wenu0000c
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  WenG_wenu0001c.root -1 1000 WenG_wenu0001c
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  WenG_wenu0000d.root -1 1000 WenG_wenu0000d
./MakeCondorFiles.csh WenG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  WenG_wenu0001d.root -1 1000 WenG_wenu0001d
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  WenG_qcd0000b.root -1 1000 WenG_qcd0000b
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  WenG_qcd0001b.root -1 1000 WenG_qcd0001b
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  WenG_qcd0002b.root -1 1000 WenG_qcd0002b
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  WenG_qcd0003b.root -1 1000 WenG_qcd0003b
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  WenG_qcd0000c.root -1 1000 WenG_qcd0000c
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  WenG_qcd0001c.root -1 1000 WenG_qcd0001c
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  WenG_qcd0000d.root -1 1000 WenG_qcd0000d
./MakeCondorFiles.csh WenG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  WenG_qcd0001d.root -1 1000 WenG_qcd0001d

nohup ./WenG_mc_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/ WenG_ZNuNuGJets.root -1 1000 >& WenG_ZNuNuGJets.txt &
nohup ./WenG_mc_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ WenG_WGJets.root -1 1000 >& WenG_WGJets.txt &
nohup ./WenG_mc_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ WenG_ZLLGJets.root -1 1000 >& WenG_ZLLGJets.txt &
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/160621_133135/0000/  WenG_TTGJets.root -1 1000 WenG_TTGJets
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/160630_154423/0000/  WenG_WZ.root -1 1000 WenG_WZ
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/160630_154436/0000/  WenG_ZZ.root -1 1000 WenG_ZZ
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_WWG/160630_154410/0000/  WenG_WWG.root -1 1000 WenG_WWG
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/160621_133030/0000/  WenG_GJets_HT-40to100.root -1 1000 WenG_GJets_HT-40to100
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/160621_133042/0000/  WenG_GJets_HT-100to200.root -1 1000 WenG_GJets_HT-100to200
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/160621_133055/0000/  WenG_GJets_HT-200to400.root -1 1000 WenG_GJets_HT-200to400
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/160621_133110/0000/  WenG_GJets_HT-400to600.root -1 1000 WenG_GJets_HT-400to600
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/160621_133123/0000/  WenG_GJets_HT-600toInf.root -1 1000 WenG_GJets_HT-600toInf
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TGJets/160823_154537/0000/  WenG_TGJets.root -1 1000 WenG_TGJets
./MakeCondorFiles.csh WenG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets_MGG-80toInf_amcatnlo/160825_000327/0000/  WenG_DiPhoton.root -1 1000 WenG_DiPhoton

./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  WmnG_data0000b.root -1 1000 WmnG_data0000b
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  WmnG_data0001b.root -1 1000 WmnG_data0001b
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  WmnG_data0002b.root -1 1000 WmnG_data0002b
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  WmnG_data0003b.root -1 1000 WmnG_data0003b
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  WmnG_data0000c.root -1 1000 WmnG_data0000c
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  WmnG_data0001c.root -1 1000 WmnG_data0001c
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  WmnG_data0000d.root -1 1000 WmnG_data0000d
./MakeCondorFiles.csh WmnG_data /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  WmnG_data0001d.root -1 1000 WmnG_data0001d
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  WmnG_wenu0000b.root -1 1000 WmnG_wenu0000b
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  WmnG_wenu0001b.root -1 1000 WmnG_wenu0001b
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  WmnG_wenu0002b.root -1 1000 WmnG_wenu0002b
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  WmnG_wenu0003b.root -1 1000 WmnG_wenu0003b
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  WmnG_wenu0000c.root -1 1000 WmnG_wenu0000c
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  WmnG_wenu0001c.root -1 1000 WmnG_wenu0001c
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  WmnG_wenu0000d.root -1 1000 WmnG_wenu0000d
./MakeCondorFiles.csh WmnG_wenu /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  WmnG_wenu0001d.root -1 1000 WmnG_wenu0001d
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0000/  WmnG_qcd0000b.root -1 1000 WmnG_qcd0000b
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0001/  WmnG_qcd0001b.root -1 1000 WmnG_qcd0001b
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0002/  WmnG_qcd0002b.root -1 1000 WmnG_qcd0002b
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016B_13TeV_2016_12p9fb_allvar_newbh_1/160720_172204/0003/  WmnG_qcd0003b.root -1 1000 WmnG_qcd0003b
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0000/  WmnG_qcd0000c.root -1 1000 WmnG_qcd0000c
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016C_13TeV_2016_12p9fb_allvar_newbh_1/160720_172250/0001/  WmnG_qcd0001c.root -1 1000 WmnG_qcd0001c
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0000/  WmnG_qcd0000d.root -1 1000 WmnG_qcd0000d
./MakeCondorFiles.csh WmnG_qcd /hdfs/store/user/gomber/SinglePhoton_2016B/SinglePhoton/crab_job_SinglePhoton2016D_13TeV_2016_12p9fb_allvar_newbh_1/160720_145740/0001/  WmnG_qcd0001d.root -1 1000 WmnG_qcd0001d

nohup ./WmnG_mc_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/ WmnG_ZNuNuGJets.root -1 1000 >& WmnG_ZNuNuGJets.txt &
nohup ./WmnG_mc_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ WmnG_WGJets.root -1 1000 >& WmnG_WGJets.txt &
nohup ./WmnG_mc_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ WmnG_ZLLGJets.root -1 1000 >& WmnG_ZLLGJets.txt &
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TTGJets_TuneCUETP8M1_13TeV-amcatnloFXFX-madspin-pythia8/crab_TTGJets/160621_133135/0000/  WmnG_TTGJets.root -1 1000 WmnG_TTGJets
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WZ_TuneCUETP8M1_13TeV-pythia8/crab_WZ/160630_154423/0000/  WmnG_WZ.root -1 1000 WmnG_WZ
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZZ_TuneCUETP8M1_13TeV-pythia8/crab_ZZ/160630_154436/0000/  WmnG_ZZ.root -1 1000 WmnG_ZZ
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WWG_TuneCUETP8M1_13TeV-amcatnlo-pythia8/crab_WWG/160630_154410/0000/  WmnG_WWG.root -1 1000 WmnG_WWG
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-40To100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-40To100/160621_133030/0000/  WmnG_GJets_HT-40to100.root -1 1000 WmnG_GJets_HT-40to100
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-100To200/160621_133042/0000/  WmnG_GJets_HT-100to200.root -1 1000 WmnG_GJets_HT-100to200
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-200To400_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-200To400/160621_133055/0000/  WmnG_GJets_HT-200to400.root -1 1000 WmnG_GJets_HT-200to400
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-400To600/160621_133110/0000/  WmnG_GJets_HT-400to600.root -1 1000 WmnG_GJets_HT-400to600
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/GJets_HT-600ToInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_GJets_HT-600ToInf/160621_133123/0000/  WmnG_GJets_HT-600toInf.root -1 1000 WmnG_GJets_HT-600toInf
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/TGJets_TuneCUETP8M1_13TeV_amcatnlo_madspin_pythia8/crab_TGJets/160823_154537/0000/  WmnG_TGJets.root -1 1000 WmnG_TGJets
./MakeCondorFiles.csh WmnG_mc /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/DiPhotonJets_MGG-80toInf_13TeV_amcatnloFXFX_pythia8/crab_DiPhotonJets_MGG-80toInf_amcatnlo/160825_000327/0000/  WmnG_DiPhoton.root -1 1000 WmnG_DiPhoton


# Consolidate output files
hadd -f ZnnG_data_all.root ZnnG_data000**.root
hadd -f ZllG_data_all.root ZllG_data000**.root
hadd -f WenG_data_all.root WenG_data000**.root
hadd -f WmnG_data_all.root WmnG_data000**.root

hadd -f ZnnG_wenu_all.root ZnnG_wenu000**.root
hadd -f ZllG_wenu_all.root ZllG_wenu000**.root
hadd -f WenG_wenu_all.root WenG_wenu000**.root
hadd -f WmnG_wenu_all.root WmnG_wenu000**.root

hadd -f ZnnG_qcd_all.root ZnnG_qcd000**.root
hadd -f WenG_qcd_all.root WenG_qcd000**.root
hadd -f WmnG_qcd_all.root WmnG_qcd000**.root

# Make fancy plots
root -l -q -b wmng_transfer_factor_plotter.C
root -l -q -b znng_transfer_factor_plotter.C
root -l -q -b znng_plotter.C
root -l -q -b zllg_plotter.C
root -l -q -b weng_plotter.C
root -l -q -b wmng_plotter.C
