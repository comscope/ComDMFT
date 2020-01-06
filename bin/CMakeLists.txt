cmake_minimum_required(VERSION 3.10.3)

project(user_interface VERSION 0.0.0 LANGUAGES NONE)

add_subdirectory(cif2matdelab)
add_subdirectory(matdelab_plot)

install(FILES cif2matdelab.py DESTINATION bin
        PERMISSIONS WORLD_READ WORLD_EXECUTE OWNER_WRITE)

install(FILES matdelab_plot.py DESTINATION bin
        PERMISSIONS WORLD_READ WORLD_EXECUTE OWNER_WRITE)
