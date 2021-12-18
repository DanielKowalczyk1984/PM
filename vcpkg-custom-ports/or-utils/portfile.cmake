
vcpkg_from_gitlab(
    GITLAB_URL https://gitlab.kuleuven.be
    OUT_SOURCE_PATH SOURCE_PATH
    REF 5e8d36decd0b374aa8c51231292711ddd85de0fb
    SHA512 6b5cabc98eb111e48eb00a168af077f6f0d60cebaf14bda4e4bc66f4e1c6d33ea34c9e078d0020cd3db607cab1a1b0058740ef1457247b4241969598983e2ed7
    REPO u0056096/or-utils
    HEAD_REF develop
)

set(ENV{GUROBI_HOME} "C:/gurobi912/win64")

vcpkg_configure_cmake(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
)

vcpkg_install_cmake()
vcpkg_fixup_cmake_targets()
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")