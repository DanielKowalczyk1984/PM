
vcpkg_from_gitlab(
    GITLAB_URL https://gitlab.kuleuven.be
    OUT_SOURCE_PATH SOURCE_PATH
    REF v0.0.4
    SHA512 1badf6530eab6fbb97348157c633979cbaa87c2c0b872b8ddf4354f58db050fa70a8a59f7f84762d0727a5006b1ffeda5d51a0218c6894c23c228f813cfb1940
    REPO u0056096/DecisionDiagramsModernCpp
)

vcpkg_cmake_configure(
    SOURCE_PATH ${SOURCE_PATH}
    PREFER_NINJA
)

vcpkg_cmake_install()
file(REMOVE_RECURSE "${CURRENT_PACKAGES_DIR}/debug")
vcpkg_copy_pdbs()