### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Add directories containing main programs
## 
######################################################################

add_subdirectory(scripts)
add_subdirectory(compare)
add_subdirectory(AnalyzeBicrystal)
add_subdirectory(AnalyzeCSL)
add_subdirectory(CreateOrientedPlane)
add_subdirectory(CreateGB)
add_subdirectory(CreateFacetGB)
add_subdirectory(SaveSteinhardtToDatabase_PerfectCrystal)
add_subdirectory(SaveSteinhardtToDatabase_Defect)
add_subdirectory(GetSteinhardt)
add_subdirectory(BondOriParam)
add_subdirectory(StructuralAnalysis_Steinhardt)
add_subdirectory(StructuralAnalysis_Steinhardt_GMM)
add_subdirectory(SaveGaussianMixture_Steinhardt)
add_subdirectory(AtomicStrain)
add_subdirectory(SearchLinkedTetrahedra)
add_subdirectory(Compute_D2Min)
add_subdirectory(AdjustStoichiometry)
add_subdirectory(IdentifyGB)
#add_subdirectory(AnalyzePolycrystal)
add_subdirectory(AnalyzeBicrystal_AndAtomicStrain)
add_subdirectory(FitAndSaveGMM)
add_subdirectory(GMMClassification)
