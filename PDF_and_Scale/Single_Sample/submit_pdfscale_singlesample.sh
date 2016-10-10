./rootcom postAnalyzer_WenG_PDFandScale WenG_pdfscale
./rootcom postAnalyzer_WenG_PDFandScale_wg WenG_pdfscale_wg
./rootcom postAnalyzer_WenG_PDFandScale_znng WenG_pdfscale_znng
./rootcom postAnalyzer_WmnG_PDFandScale WmnG_pdfscale
./rootcom postAnalyzer_WmnG_PDFandScale_wg WmnG_pdfscale_wg
./rootcom postAnalyzer_WmnG_PDFandScale_znng WmnG_pdfscale_znng
./rootcom postAnalyzer_ZllG_PDFandScale ZllG_pdfscale
./rootcom postAnalyzer_ZllG_PDFandScale_wg ZllG_pdfscale_wg
./rootcom postAnalyzer_ZllG_PDFandScale_znng ZllG_pdfscale_znng
./rootcom postAnalyzer_ZnnG_PDFandScale ZnnG_pdfscale
./rootcom postAnalyzer_ZnnG_PDFandScale_wg ZnnG_pdfscale_wg
./rootcom postAnalyzer_ZnnG_PDFandScale_znng ZnnG_pdfscale_znng

nohup ./ZllG_pdfscale_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ ZllG_pdfscale_ZLLGJets.root -1 1000 >& ZllG_pdfscale_ZLLGJets.txt &
nohup ./ZllG_pdfscale_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ ZllG_pdfscale_WGJets.root -1 1000 >& ZllG_pdfscale_WGJets.txt &

nohup ./WenG_pdfscale_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ WenG_pdfscale_WGJets.root -1 1000 >& WenG_pdfscale_WGJets.txt &
nohup ./WenG_pdfscale_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ WenG_pdfscale_ZLLGJets.root -1 1000 >& WenG_pdfscale_ZLLGJets.txt &
nohup ./WmnG_pdfscale_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ WmnG_pdfscale_WGJets.root -1 1000 >& WmnG_pdfscale_WGJets.txt &
nohup ./WmnG_pdfscale_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ WmnG_pdfscale_ZLLGJets.root -1 1000 >& WmnG_pdfscale_ZLLGJets.txt &

nohup ./ZnnG_pdfscale_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZNuNuGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZNuNuGJets/160621_132818/0000/ ZnnG_pdfscale_ZNuNuGJets.root -1 1000 >& ZnnG_pdfscale_ZNuNuGJets.txt &
nohup ./ZnnG_pdfscale_znng /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/ZLLGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_ZLLGJets/160621_132831/0000/ ZnnG_pdfscale_ZLLGJets.root -1 1000 >& ZnnG_pdfscale_ZLLGJets.txt &
nohup ./ZnnG_pdfscale_wg /hdfs/store/user/jjbuch/ntuples80X/MiniAODv2_CustomVariables/WGJets_MonoPhoton_PtG-130_TuneCUETP8M1_13TeV-madgraph/crab_WGJets/160621_132844/0000/ ZnnG_pdfscale_WGJets.root -1 1000 >& ZnnG_pdfscale_WGJets.txt &
