// file: sdl_gl_sphere.c
#include <SDL.h>
#include <SDL_opengl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char** argv) {
    if (SDL_Init(SDL_INIT_VIDEO) != 0) {
        fprintf(stderr, "SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    // OpenGL özellikleri (eski tip, compatibility profile)
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
    SDL_Window* win = SDL_CreateWindow("SDL2 + OpenGL Sphere",
                                       SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED,
                                       800, 600, SDL_WINDOW_OPENGL | SDL_WINDOW_SHOWN);
    if (!win) {
        fprintf(stderr, "SDL_CreateWindow error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_GLContext ctx = SDL_GL_CreateContext(win);
    if (!ctx) {
        fprintf(stderr, "SDL_GL_CreateContext error: %s\n", SDL_GetError());
        SDL_DestroyWindow(win);
        SDL_Quit();
        return 1;
    }

    // Basit GL ayarları
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat lightPos[] = { 1.0f, 1.0f, 2.0f, 0.0f };
    glLightfv(GL_LIGHT0, GL_POSITION, lightPos);

    int running = 1;
    SDL_Event ev;
    float angle = 0.0f;

    while (running) {
        while (SDL_PollEvent(&ev)) {
            if (ev.type == SDL_QUIT) running = 0;
            if (ev.type == SDL_KEYDOWN && ev.key.keysym.sym == SDLK_ESCAPE) running = 0;
        }

        int w, h;
        SDL_GetWindowSize(win, &w, &h);
        glViewport(0, 0, w, h);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, (double)w / (double)h, 0.1, 100.0);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        gluLookAt(0.0, 0.0, 4.0,  // kamera pos
                  0.0, 0.0, 0.0,  // hedef
                  0.0, 1.0, 0.0); // up

        glClearColor(0.08f, 0.08f, 0.12f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // döndürme ile animasyon
        glPushMatrix();
        glRotatef(angle, 0.3f, 1.0f, 0.0f);

        // Materyal
        GLfloat mat_ambient[] = {0.2f, 0.2f, 0.6f, 1.0f};
        GLfloat mat_diffuse[] = {0.2f, 0.6f, 0.9f, 1.0f};
        GLfloat mat_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
        GLfloat mat_shininess[] = {50.0f};
        glMaterialfv(GL_FRONT, GL_AMBIENT, mat_ambient);
        glMaterialfv(GL_FRONT, GL_DIFFUSE, mat_diffuse);
        glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
        glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);

        GLUquadric* q = gluNewQuadric();
        gluQuadricNormals(q, GLU_SMOOTH);
        gluSphere(q, 1.0, 40, 40); // radius=1.0, slices=40, stacks=40
        gluDeleteQuadric(q);

        glPopMatrix();

        SDL_GL_SwapWindow(win);

        angle += 0.6f; // döndürme hızı
        SDL_Delay(16); // ~60 FPS
    }

    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
