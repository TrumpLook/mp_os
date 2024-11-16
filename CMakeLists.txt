cmake_minimum_required(VERSION 3.26)
project(cw_os)

set(CMAKE_CXX_STANDARD 17)

add_executable(cw_os
        main.cpp
        C:/OpSys_cw/cw/database/src/data_base.cpp
        C:/OpSys_cw/cw/schema/src/schema.cpp
        C:/OpSys_cw/cw/schemas_pool/src/schemas_pool.cpp
        C:/OpSys_cw/cw/table/src/table.cpp
        C:/OpSys_cw/b-tree/b_tree.h
        C:/OpSys_cw/cw/user_data/src/user_data.cpp
        C:/OpSys_cw/allocator/allocator_boundary_tags/src/allocator_boundary_tags.cpp
        C:/OpSys_cw/main.cpp
        C:/OpSys_cw/allocator/allocator/src/allocator_guardant.cpp
        C:/OpSys_cw/logger/client_logger/src/client_logger.cpp
        C:/OpSys_cw/logger/client_logger/src/client_logger_builder.cpp
        C:/OpSys_cw/logger/logger/src/logger.cpp
        C:/OpSys_cw/logger/logger/src/logger_builder.cpp
        C:/OpSys_cw/logger/logger/src/logger_guardant.cpp
        logger/client_logger/src/client_logger_builder.cpp
        C:/OpSys_cw/cw/string_pool/src/string_pool.cpp
)

