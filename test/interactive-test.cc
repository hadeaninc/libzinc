#include <GLFW/glfw3.h>
#define GLFW_EXPOSE_NATIVE_X11
#define GLFW_EXPOSE_NATIVE_GLX
#include <GLFW/glfw3native.h>

#include <cairo/cairo.h>
#include <cairo/cairo-gl.h>

#include <libzinc/zinc.hh>

int main() {
    glfwInit();

    uint32_t width = 800, height = 800;

    GLFWwindow* window = glfwCreateWindow(width, height, "interactive test", NULL, NULL);
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    Display* x11_display = glfwGetX11Display();
    GLXContext glx_context = glfwGetGLXContext(window);
    Window x11_window = glfwGetX11Window(window);

    cairo_device_t* device = cairo_glx_device_create(x11_display, glx_context);
    cairo_surface_t* surface = cairo_gl_surface_create_for_window(device, x11_window, width, height);
    cairo_device_destroy(device);
    cairo_t* cr = cairo_create(surface);

    while (!glfwWindowShouldClose(window)) {
        double mouse_x, mouse_y;
        glfwGetCursorPos(window, &mouse_x, &mouse_y);

        cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
        cairo_paint(cr);

        cairo_save(cr);

        float scale = 16.0;
        cairo_translate(cr, scale, scale);
        cairo_scale(cr, scale, scale);

        {
            cairo_move_to(cr, 0, 0);
            for (uint64_t m = 1; m < 1024; m++) {
                auto p = morton_code<2, 32>::decode({m});
                uint32_t x = p[0];
                uint32_t y = p[1];
                cairo_line_to(cr, x, y);
            }
        }

        cairo_restore(cr);

        cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);
        cairo_set_line_width(cr, 1.0);
        cairo_stroke(cr);

        cairo_gl_surface_swapbuffers(surface);
        glfwSwapBuffers(window);

        glfwPollEvents();
        if (glfwGetKey(window, GLFW_KEY_Q) == GLFW_PRESS) {
            glfwSetWindowShouldClose(window, 1);
        }
    }

    cairo_surface_destroy(surface);
    cairo_destroy(cr);

    glfwTerminate();
}
