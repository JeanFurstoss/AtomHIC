### CMakeLists.txt --- 
##
######################################################################
## 
### Commentary: Add directories containing main programs
## 
######################################################################

add_subdirectory(CreateCrystal)
add_subdirectory(AnalyzeBicrystal)
add_subdirectory(AnalyzeCSL)
add_subdirectory(AnalyzeBicrystal_ARGS)
add_subdirectory(OrderParam)
add_subdirectory(CreateOrientedPlane)
add_subdirectory(CreateGB)
add_subdirectory(CreateFacetGB)
add_subdirectory(ComputeStrainTensor_FromNum)
add_subdirectory(SaveSteinhardtToDatabase_PerfectCrystal)
add_subdirectory(SaveSteinhardtToDatabase_Defect)
add_subdirectory(OrderParam_Steinhardt)
add_subdirectory(ComputeSteinhardtParams)
add_subdirectory(StructuralAnalysis_Steinhardt)

