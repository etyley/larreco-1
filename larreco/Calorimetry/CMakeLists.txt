art_make(LIB_LIBRARIES
           lardataobj_RecoBase
           lardata_ArtDataHelper
           larcorealg_Geometry
           lardataobj_AnalysisBase
           ${ART_FRAMEWORK_SERVICES_REGISTRY}
           ${ART_UTILITIES}
           canvas
           ${MF_MESSAGELOGGER}
           ${FHICLCPP}
           cetlib_except
           ROOT::Core
           ROOT::Physics
         MODULE_LIBRARIES
           larcorealg_Geometry
           larreco_Calorimetry
           ROOT::Core
           ROOT::Hist
           ROOT::Physics
           ${MF_MESSAGELOGGER}
         )

install_headers()
install_fhicl()
install_source()
