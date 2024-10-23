add_executable(ConInter_demo)
add_executable(ConInter::demo ALIAS ConInter_demo)
target_sources(ConInter_demo PRIVATE interpolation_demo.cxx)
target_link_libraries(ConInter_demo PRIVATE ConInter)

add_executable(ConInter_demo_largestrain)
add_executable(ConInter::demo_largestrain ALIAS ConInter_demo_largestrain)
target_sources(ConInter_demo_largestrain
               PRIVATE interpolation_demo_largestrain.cxx)
target_link_libraries(ConInter_demo_largestrain PRIVATE ConInter)
