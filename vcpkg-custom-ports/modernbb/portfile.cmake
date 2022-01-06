vcpkg_from_gitlab(
  GITLAB_URL
  https://gitlab.kuleuven.be
  OUT_SOURCE_PATH
  SOURCE_PATH
  REF
  v0.1.2
  SHA512
  9cc87c073ad3e3b337a14327aae1a2aeffda1e40c3102c27e91d291f1e1124aa4fa71f428d5cf4f48ceb7b668a224891de2b8571bf05a04b6a905fdd450fa4f3
  REPO
  u0056096/branch-and-bound
)

vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH} PREFER_NINJA)

vcpkg_cmake_install()
vcpkg_cmake_config_fixup(PACKAGE_NAME "branch-and-bound")
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_copy_pdbs()
