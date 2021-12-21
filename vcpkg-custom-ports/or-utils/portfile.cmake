vcpkg_from_gitlab(
  GITLAB_URL
  https://gitlab.kuleuven.be
  OUT_SOURCE_PATH
  SOURCE_PATH
  REF
  v1.0.0
  SHA512
  19dc2bfa7e291a90fe26ff7e9afb324457997adccabb389dace2a2f49f5a2beb3c854aa9cefcbfa48eab7e351417885bb430822a1dc12cae910ccb4790f3aca1
  REPO
  u0056096/or-utils
  HEAD_REF
  develop
)

if(VCPKG_TARGET_IS_WINDOWS)
  set(ENV{GUROBI_HOME} "C:/gurobi912/win64")
elseif(VCPKG_IS_LINUX)
  set(ENV{GUROBI_HOME} "/opt/gurobi/linux64")
endif()
vcpkg_cmake_configure(SOURCE_PATH ${SOURCE_PATH} PREFER_NINJA)

vcpkg_cmake_install()
file(
  INSTALL "${SOURCE_PATH}/LICENSE.md"
  DESTINATION "${CURRENT_PACKAGES_DIR}/share/${PORT}"
  RENAME copyright
)
vcpkg_cmake_config_fixup()
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug/include")
vcpkg_copy_pdbs()
