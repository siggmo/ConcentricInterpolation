# The source files are grouped corresponding to the Makefile targets in the
# original project

# ##############################################################################
# lib (src/Makefile:23), consisting of the following targets:
# ##############################################################################

# CI (src/Makefile:34)
target_sources(ConInter PRIVATE concentric_interpolation.cxx)

# RI (src/Makefile:37)
target_sources(ConInter PRIVATE radial_interpolation.cxx)

# TI (src/Makefile:40)
target_sources(ConInter PRIVATE tangential_interpolation.cxx)

# data (src/Makefile:43)
target_sources(ConInter PRIVATE data.cxx)

# util (src/Makefile:46)
target_sources(ConInter PRIVATE util.cxx)
