load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

http_archive(
    name = "platforms",
    urls = [
        "https://mirror.bazel.build/github.com/bazelbuild/platforms/releases/download/0.0.6/platforms-0.0.6.tar.gz",
        "https://github.com/bazelbuild/platforms/releases/download/0.0.6/platforms-0.0.6.tar.gz",
    ],
    sha256 = "5308fc1d8865406a49427ba24a9ab53087f17f5266a7aabbfc28823f3916e1ca",
)

http_archive(
    name = "com_google_googletest",
    sha256 = "ab78fa3f912d44d38b785ec011a25f26512aaedc5291f51f3807c592b506d33a",
    strip_prefix = "googletest-58d77fa8070e8cec2dc1ed015d66b454c8d78850",
    url = "https://github.com/google/googletest/archive/58d77fa8070e8cec2dc1ed015d66b454c8d78850.zip",
)

# Required for testing compatibility with TF Quantum:
# https://github.com/tensorflow/quantum
http_archive(
    name = "org_tensorflow",
    sha256 = "e52cda3bae45f0ae0fccd4055e9fa29892b414f70e2df94df9a3a10319c75fff",
    strip_prefix = "tensorflow-2.11.0",
    urls = [
        "https://github.com/tensorflow/tensorflow/archive/refs/tags/v2.11.0.zip",
    ],
)

load("@org_tensorflow//tensorflow:workspace3.bzl", "workspace")

workspace()

load("@org_tensorflow//tensorflow:workspace2.bzl", "workspace")

workspace()

load("@org_tensorflow//tensorflow:workspace1.bzl", "workspace")

workspace()

load("@org_tensorflow//tensorflow:workspace0.bzl", "workspace")

workspace()


EIGEN_COMMIT = "3bb6a48d8c171cf20b5f8e48bfb4e424fbd4f79e"
EIGEN_SHA256 = "eca9847b3fe6249e0234a342b78f73feec07d29f534e914ba5f920f3e09383a3"


http_archive(
    name = "eigen",
    build_file_content = """
cc_library(
  name = "eigen3",
  textual_hdrs = glob(["Eigen/**", "unsupported/**"]),
  visibility = ["//visibility:public"],
)
    """,
    sha256 = EIGEN_SHA256,
        strip_prefix = "eigen-{commit}".format(commit = EIGEN_COMMIT),
        urls = [
            "https://storage.googleapis.com/mirror.tensorflow.org/gitlab.com/libeigen/eigen/-/archive/{commit}/eigen-{commit}.tar.gz".format(commit = EIGEN_COMMIT),
            "https://gitlab.com/libeigen/eigen/-/archive/{commit}/eigen-{commit}.tar.gz".format(commit = EIGEN_COMMIT),
        ],
)

load("//third_party/cuquantum:cuquantum_configure.bzl", "cuquantum_configure")

cuquantum_configure(name = "local_config_cuquantum")

# External dependency: xtl (dependency for xtensor)
http_archive(
    name = "xtl",
    urls = ["https://github.com/xtensor-stack/xtl/archive/refs/heads/master.tar.gz"],
    strip_prefix = "xtl-master",
    build_file_content = """
cc_library(
    name = "xtl",
    hdrs = glob(["xtl/**/*.hpp"]),
    includes = ["xtl"],
    visibility = ["//visibility:public"],
)
    """,
)

# External dependency: xtensor
http_archive(
    name = "xtensor",
    urls = ["https://github.com/xtensor-stack/xtensor/archive/refs/heads/master.tar.gz"],
    strip_prefix = "xtensor-master",
    build_file_content = """
cc_library(
    name = "xtensor",
    hdrs = glob(["xtensor/**/*.hpp"]),
    includes = ["xtensor"],
    deps = [":xtl"],
    visibility = ["//visibility:public"],
)
    """,
)

# External dependency: xtensor-blas
http_archive(
    name = "xtensor-blas",
    urls = ["https://github.com/xtensor-stack/xtensor-blas/archive/refs/heads/master.tar.gz"],
    strip_prefix = "xtensor-blas-master",
    build_file_content = """
cc_library(
    name = "xtensor-blas",
    hdrs = glob(["xtensor-blas/**/*.hpp"]),
    includes = ["xtensor-blas"],
    deps = [":xtensor"],
    visibility = ["//visibility:public"],
)
    """,
)
