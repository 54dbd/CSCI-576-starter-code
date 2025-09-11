-- wxWidgets 子模块 xmake.lua 示例
-- 这里只做最基础的头文件和库目录导出，实际项目建议用 CMake 构建 wxWidgets

-- 导出头文件目录，供父项目 add_includedirs
add_includedirs("include")
add_includedirs("include/msw")

-- 可选: 导出库目录，供父项目 add_linkdirs
-- add_linkdirs("lib/vc_x64_dll")

-- 可选: 定义 wxWidgets 相关 target（如有源码/自定义模块）
-- target("wxWidgets")
--     set_kind("static")
--     add_files("src/*.cpp")
--     add_includedirs("include")
--     add_includedirs("include/msw")
