cmake_minimum_required(VERSION 3.10.3)

project(matdelab_plot VERSION 0.0.0 LANGUAGES NONE)

install(FILES band_structure.py DESTINATION bin/matdelab_plot)

install(FILES density_of_states.py DESTINATION bin/matdelab_plot)

install(FILES figure_bd.py DESTINATION bin/matdelab_plot)
