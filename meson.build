project('zinc', 'cpp',
  default_options: ['b_sanitize=address', 'b_lundef=false']
)

add_project_arguments(
  '-mbmi2',
  '-Werror', '-Wall', '-Wextra', '-Wpedantic',
  '-Wno-unused-parameter',
  '-fsanitize=implicit-conversion', '-fsanitize=integer', '-fsanitize=undefined',
  '-std=c++17',
  language: 'cpp',
)

cxx = meson.get_compiler('cpp')
m_dep = cxx.find_library('m', required : false)

incdir = include_directories('libzinc')

glfw_dep = dependency('glfw3')
cairo_dep = dependency('cairo')
interactive_test = executable('interactive-test', 'test/interactive-test.cc', include_directories: incdir, dependencies: [m_dep, glfw_dep, cairo_dep])

zinc_test = executable('zinc-test', 'test/zinc-test.cc', include_directories: incdir, dependencies: m_dep)
test('zinc-test', zinc_test)

install_subdir('libzinc', install_dir: 'include')
