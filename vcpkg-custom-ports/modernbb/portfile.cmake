
vcpkg_from_gitlab(
    GITLAB_URL https://gitlab.kuleuven.be
    OUT_SOURCE_PATH SOURCE_PATH
    REF v0.1.1
    SHA512 11b7b6ec26045ffbe0c5c84e79ddee3baed8183829a3cc452fe590f7270afeef703c319ca8d7528f96ec3f8e6436586bf82c25abf70c3018ef959d1cccaf9a41
    REPO u0056096/branch-and-bound
)

vcpkg_cmake_configure(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME "branch-and-bound")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_copy_pdbs()