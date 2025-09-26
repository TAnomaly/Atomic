// Professional Atomic Simulation with Real Orbitals
#include <SDL.h>
#include <SDL_opengl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Configuration constants
#define NUM_ELECTRONS 4
#define TRAIL_LEN 150
#define NUCLEUS_RADIUS 0.8f

// Visual quality settings for perfect spheres
#define SPHERE_QUALITY 48  // Higher quality for smoother spheres
#define GLOW_INTENSITY 0.3f
#define TRAIL_ALPHA_START 0.9f
#define TRAIL_ALPHA_END 0.0f
#define TRAIL_FADE_SPEED 0.015f

typedef struct {
    float x, y, z;
} Vec3;

typedef struct {
    Vec3 position;
    Vec3 color;
    float alpha;
} TrailPoint;

typedef struct {
    float radius;
    float speed;
    float angle_x, angle_y, angle_z;
    Vec3 color;
} Electron;

// Global state
TrailPoint electronTrail[NUM_ELECTRONS][TRAIL_LEN];
int trailIndex[NUM_ELECTRONS] = {0};
Electron electrons[NUM_ELECTRONS];
float global_time = 0.0f;

// Initialize simple atom structure
void initializeAtom() {
    srand(time(NULL));
    
    // Initialize electrons with different orbital radii and colors
    for(int i = 0; i < NUM_ELECTRONS; i++) {
        electrons[i].radius = 2.5f + i * 1.2f; // Different orbital radii
        electrons[i].speed = 0.8f + ((float)rand()/RAND_MAX) * 0.6f; // Random speeds
        electrons[i].angle_x = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
        electrons[i].angle_y = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
        electrons[i].angle_z = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
        
        // Different colors for each electron
        switch(i) {
            case 0: electrons[i].color = (Vec3){0.2f, 0.8f, 1.0f}; break; // Cyan
            case 1: electrons[i].color = (Vec3){1.0f, 0.3f, 0.3f}; break; // Red
            case 2: electrons[i].color = (Vec3){0.3f, 1.0f, 0.3f}; break; // Green
            case 3: electrons[i].color = (Vec3){1.0f, 0.8f, 0.2f}; break; // Yellow
        }
    }
    
    // Initialize all trail points
    for(int i = 0; i < NUM_ELECTRONS; i++) {
        for(int j = 0; j < TRAIL_LEN; j++) {
            electronTrail[i][j].position = (Vec3){0, 0, 0};
            electronTrail[i][j].color = electrons[i].color;
            electronTrail[i][j].alpha = 0.0f;
        }
    }
}

// Perfect 3D sphere rendering with smooth surface
void drawPerfectSphere(float radius, GLfloat r, GLfloat g, GLfloat b, float alpha) {
    // Set up material properties for realistic sphere
    GLfloat mat_ambient[]  = { r*0.2f, g*0.2f, b*0.2f, alpha };
    GLfloat mat_diffuse[]  = { r*0.7f, g*0.7f, b*0.7f, alpha };
    GLfloat mat_specular[] = { 0.9f, 0.9f, 0.9f, alpha };
    GLfloat mat_shininess  = 100.0f;
    GLfloat mat_emission[] = { r*0.05f, g*0.05f, b*0.05f, alpha };
    
    glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, mat_ambient);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, mat_diffuse);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, mat_specular);
    glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, mat_shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_EMISSION, mat_emission);

    // Enable smooth shading
    glShadeModel(GL_SMOOTH);
    
    // Create high-quality sphere
    GLUquadric* sphere = gluNewQuadric();
    gluQuadricNormals(sphere, GLU_SMOOTH);
    gluQuadricTexture(sphere, GL_FALSE);
    gluQuadricDrawStyle(sphere, GLU_FILL);
    
    // Draw perfect sphere
    gluSphere(sphere, radius, SPHERE_QUALITY, SPHERE_QUALITY);
    gluDeleteQuadric(sphere);
}

// Draw soft glowing effect around spheres
void drawSoftGlow(float radius, GLfloat r, GLfloat g, GLfloat b) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING);
    glDepthMask(GL_FALSE); // Don't write to depth buffer for glow
    
    // Multi-layer glow effect
    for(int i = 0; i < 4; i++) {
        float glow_radius = radius * (1.3f + i * 0.2f);
        float alpha = GLOW_INTENSITY * (1.0f - (float)i / 4.0f);
        
        glColor4f(r, g, b, alpha);
        
        GLUquadric* glow_sphere = gluNewQuadric();
        gluQuadricNormals(glow_sphere, GLU_SMOOTH);
        gluQuadricDrawStyle(glow_sphere, GLU_FILL);
        gluSphere(glow_sphere, glow_radius, 16, 16);
        gluDeleteQuadric(glow_sphere);
    }
    
    glDepthMask(GL_TRUE);
    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
}

// Perfect nucleus rendering with proper 3D sphere
void drawNucleus() {
    glPushMatrix();
    
    // Add subtle pulsing effect for energy
    float pulse = 1.0f + 0.03f * sinf(global_time * 1.5f);
    glScalef(pulse, pulse, pulse);
    
    // Draw glowing aura first
    drawSoftGlow(NUCLEUS_RADIUS, 1.0f, 0.4f, 0.1f);
    
    // Draw perfect sphere nucleus
    drawPerfectSphere(NUCLEUS_RADIUS, 0.95f, 0.6f, 0.2f, 1.0f);
    
    glPopMatrix();
}

// Enhanced trail rendering with dynamic fade
void drawElectronTrail(int electron_id) {
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_LIGHTING);
    glLineWidth(3.0f);
    
    glBegin(GL_LINE_STRIP);
    for(int i = 0; i < TRAIL_LEN; i++) {
        TrailPoint* point = &electronTrail[electron_id][i];
        if(point->alpha > 0.0f) {
            glColor4f(point->color.x, point->color.y, point->color.z, point->alpha);
            glVertex3f(point->position.x, point->position.y, point->position.z);
        }
    }
    glEnd();
    
    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
}

// Update trail system
void updateTrail(int electron_id, Vec3 position) {
    // Add new trail point
    TrailPoint* newPoint = &electronTrail[electron_id][trailIndex[electron_id]];
    newPoint->position = position;
    newPoint->color = electrons[electron_id].color;
    newPoint->alpha = TRAIL_ALPHA_START;
    
    trailIndex[electron_id] = (trailIndex[electron_id] + 1) % TRAIL_LEN;
    
    // Fade all existing trail points
    for(int i = 0; i < TRAIL_LEN; i++) {
        electronTrail[electron_id][i].alpha -= TRAIL_FADE_SPEED;
        if(electronTrail[electron_id][i].alpha < 0.0f) {
            electronTrail[electron_id][i].alpha = 0.0f;
        }
    }
}

// Calculate 3D orbital positions with smooth movement
Vec3 calculateElectronPosition(int electron_id, float time) {
    Electron* e = &electrons[electron_id];
    Vec3 pos = {0};
    
    // Update angles for 3D movement
    e->angle_x += time * e->speed * 0.7f;
    e->angle_y += time * e->speed * 0.5f;
    e->angle_z += time * e->speed * 0.3f;
    
    // Create 3D orbital motion with Lissajous curves for interesting patterns
    pos.x = e->radius * cosf(e->angle_x) * cosf(e->angle_z * 0.3f);
    pos.y = e->radius * sinf(e->angle_y) * 0.8f;
    pos.z = e->radius * sinf(e->angle_x) * sinf(e->angle_z * 0.5f);
    
    // Add small random variation for more natural look
    float variation = 0.05f;
    pos.x += (((float)rand()/RAND_MAX) - 0.5f) * variation;
    pos.y += (((float)rand()/RAND_MAX) - 0.5f) * variation;
    pos.z += (((float)rand()/RAND_MAX) - 0.5f) * variation;
    
    return pos;
}

int main(int argc, char** argv) {
    if(SDL_Init(SDL_INIT_VIDEO)!=0){fprintf(stderr,"SDL_Init error: %s\n",SDL_GetError());return 1;}
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION,2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION,1);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLEBUFFERS, 1);
    SDL_GL_SetAttribute(SDL_GL_MULTISAMPLESAMPLES, 4); // 4x MSAA

    SDL_Window* win = SDL_CreateWindow("Professional Atomic Simulation - Carbon Atom",
        SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 1200, 900, SDL_WINDOW_OPENGL|SDL_WINDOW_SHOWN|SDL_WINDOW_RESIZABLE);
    if(!win){fprintf(stderr,"SDL_CreateWindow error: %s\n",SDL_GetError());SDL_Quit();return 1;}
    SDL_GLContext ctx = SDL_GL_CreateContext(win);
    if(!ctx){fprintf(stderr,"SDL_GL_CreateContext error: %s\n",SDL_GetError());SDL_DestroyWindow(win);SDL_Quit();return 1;}

    // Initialize atom structure
    initializeAtom();

    // Enhanced OpenGL setup for perfect spheres
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glEnable(GL_LIGHTING);
    glEnable(GL_NORMALIZE);
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    
    // Enable smooth shading globally
    glShadeModel(GL_SMOOTH);
    
    // Enable back face culling for better performance
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    
    // Improve sphere quality
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
    glHint(GL_POLYGON_SMOOTH_HINT, GL_NICEST);
    
    // Multi-light setup for better visual quality
    glEnable(GL_LIGHT0); // Main light
    GLfloat light0_pos[] = {3.0f, 4.0f, 5.0f, 1.0f};
    GLfloat light0_ambient[] = {0.2f, 0.2f, 0.3f, 1.0f};
    GLfloat light0_diffuse[] = {0.8f, 0.8f, 1.0f, 1.0f};
    GLfloat light0_specular[] = {1.0f, 1.0f, 1.0f, 1.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, light0_pos);
    glLightfv(GL_LIGHT0, GL_AMBIENT, light0_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light0_specular);
    
    glEnable(GL_LIGHT1); // Fill light
    GLfloat light1_pos[] = {-2.0f, -1.0f, 3.0f, 1.0f};
    GLfloat light1_ambient[] = {0.1f, 0.1f, 0.1f, 1.0f};
    GLfloat light1_diffuse[] = {0.4f, 0.4f, 0.6f, 1.0f};
    glLightfv(GL_LIGHT1, GL_POSITION, light1_pos);
    glLightfv(GL_LIGHT1, GL_AMBIENT, light1_ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
    
    // Global ambient light
    GLfloat global_ambient[] = {0.1f, 0.1f, 0.15f, 1.0f};
    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, global_ambient);

    int running = 1;
    SDL_Event ev;
    float camera_angle = 0.0f;
    float camera_distance = 12.0f;
    
    printf("Controls:\n");
    printf("ESC - Exit\n");
    printf("Mouse Wheel - Zoom in/out\n");
    printf("Arrow Keys - Rotate view\n");
    printf("Space - Reset view\n\n");

    while(running) {
        while(SDL_PollEvent(&ev)) {
            if(ev.type == SDL_QUIT) running = 0;
            if(ev.type == SDL_KEYDOWN) {
                switch(ev.key.keysym.sym) {
                    case SDLK_ESCAPE: running = 0; break;
                    case SDLK_LEFT: camera_angle -= 0.1f; break;
                    case SDLK_RIGHT: camera_angle += 0.1f; break;
                    case SDLK_UP: camera_distance -= 0.5f; break;
                    case SDLK_DOWN: camera_distance += 0.5f; break;
                    case SDLK_SPACE: camera_angle = 0.0f; camera_distance = 12.0f; break;
                }
                if(camera_distance < 5.0f) camera_distance = 5.0f;
                if(camera_distance > 25.0f) camera_distance = 25.0f;
            }
        }

        global_time += 0.016f; // ~60 FPS timing

        int w, h;
        SDL_GetWindowSize(win, &w, &h);
        glViewport(0, 0, w, h);
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        gluPerspective(45.0, (double)w/(double)h, 0.1, 100.0);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        
        // Dynamic camera movement
        float cam_x = camera_distance * sinf(camera_angle);
        float cam_z = camera_distance * cosf(camera_angle);
        gluLookAt(cam_x, 3.0, cam_z, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        
        // Dark space background with subtle gradient
        glClearColor(0.02f, 0.02f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Render nucleus with proper structure
        drawNucleus();

        // Render electrons with 3D motion and trails
        for(int i = 0; i < NUM_ELECTRONS; i++) {
            // Calculate new position
            Vec3 pos = calculateElectronPosition(i, 0.016f);
            
            // Update trail system
            updateTrail(i, pos);
            
            // Draw fading trail
            drawElectronTrail(i);

            // Draw perfect sphere electron with soft glow
            glPushMatrix();
            glTranslatef(pos.x, pos.y, pos.z);
            Vec3 color = electrons[i].color;
            
            // Draw soft glow around electron
            drawSoftGlow(0.08f, color.x, color.y, color.z);
            
            // Draw perfect electron sphere
            drawPerfectSphere(0.1f, color.x, color.y, color.z, 0.95f);
            
            glPopMatrix();
        }

        SDL_GL_SwapWindow(win);
        SDL_Delay(16);
    }

    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
