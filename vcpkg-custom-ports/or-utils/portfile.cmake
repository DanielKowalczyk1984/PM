
vcpkg_from_gitlab(
    GITLAB_URL https://gitlab.kuleuven.be
    OUT_SOURCE_PATH SOURCE_PATH
    REF 22b44eb1155428973f57937a989692ad55ba7ef5
    SHA512 310ce569ff47b270bede1181ed5d4f445cff1f47b3121b1e6c4a97c8684c358dafae10d1d21e3611041adc35e1123b59f71e21cc2b3e74168e36d7fed82c8a99
    REPO u0056096/or-utils
    HEAD_REF develop
)

set(ENV{GUROBI_HOME} "C:/gurobi912/win64")
# if (VCPKG_LVIBRARY_LINKAGE STREQUAL dynamic)
#     message(STATUS "Warning: Dynamic building not supported yet. Building static.")
#     set(VCPKG_LIBRARY_LINKAGE static)
# endif()
vcpkg_cmake_configure(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup()
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_copy_pdbs()