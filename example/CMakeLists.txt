add_executable(gaussquad_example)
zero_get_files(tmp .)
target_sources(gaussquad_example PRIVATE ${tmp})
target_link_libraries(gaussquad_example PRIVATE gaussquad)

zero_check_target(gaussquad_example)
