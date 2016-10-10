**Control_Region_Analyzers**: Compile and run commands are listed in submit.sh, along with the hadd commands to merge the output root files at the end. CMSSW_8_0_8 with ewk_corr.root, MakeCondorFiles.csh, and rootcom are required.

**PDF_and_Scale**: There are two separate cases to consider, each having their own corresponding folder: systematics on transfer factors, and systematics on the expected yield from just a single MC sample. MC samples which could potentially be considered include ZNuNuGJets, WGJets, and ZLLGJets (the only ones with PDF weights included in the appropriate format).

**PDF_and_Scale/Single_Sample**: Compile and run commands are listed in submit_pdfscale_singlesample.sh. CMSSW_8_0_8 with ewk_corr.root and rootcom are required. Get the final "data card ready" histograms by running background_systematics_plotter.C.

**PDF_and_Scale/Transfer_Factor**: Not yet implemented.

In all instances, the suffix "\_wg" on an analyzer means that EWK and NNLO corrections appropriate to WGJets have been applied, and similarly for "\_znng". Corrections particular to ZLLGJets have not yet been supplied; a _temporary_ workaround is to use ZnnG corrections for this sample.
