# 最低版本
cmake_minimum_required (VERSION 2.8)

# 项目名称
project (OCM)
set(CMAKE_BUILD_TYPE "Release")
set(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O2 -static  -g")
set(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O2 -static")

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_C_STANDARD 11)
# 引入头文件路径
set(source_dir ${PROJECT_SOURCE_DIR}/../../)
include_directories(${source_dir} ${source_dir}/GOAL ${source_dir}/GOAL/Common ${source_dir}/GOAL/Optimization/Interface ${source_dir}/GOAL/Optimization/OCM ${source_dir}/GOAL/System ${source_dir}/GOAL/ThirdParty/CgDescent)

# 获取某目录下的所有源文件
aux_source_directory(${source_dir}/GOAL files)
aux_source_directory(${source_dir}/GOAL/Optimization/OCM files1)
aux_source_directory(${source_dir}/Test files2)

# 生产可执行文件
add_executable (OCM ${files} ${files1} ${files2}  )
