file(GLOB_RECURSE CUDA_CUBBYFLOW_SOURCE *.h *.hpp *.cpp *.cu)
add_library(CUDA_CubbyFlow STATIC ${CUDA_CUBBYFLOW_SOURCE})
target_link_libraries(CUDA_CubbyFlow PUBLIC ${CUDA_LIBRARIES} tbb)
target_link_directories(
        CUDA_CubbyFlow
        PUBLIC
        ${Houdini_PATH}/custom/houdini/dsolib
)
target_include_directories(
        CUDA_CubbyFlow
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/tinyobj
        ${CMAKE_CURRENT_SOURCE_DIR}/flatbuffers
        ${CMAKE_CURRENT_SOURCE_DIR}/Core
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/../
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
        ${_houdini_include_dir}
)
target_compile_definitions(
        CUDA_CubbyFlow
        PUBLIC
        CUBBYFLOW_TASKING_TBB
)
target_compile_definitions(
        CUDA_CubbyFlow
        PUBLIC
        CUBBYFLOW_USE_CUDA
)