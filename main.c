// Professional Atomic Simulation with Real Orbitals
#include <SDL.h>
#include <SDL_opengl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>

// Configuration constants
#define NUM_MOLECULES 4
#define NUM_ELECTRONS_PER_MOLECULE 3
#define TOTAL_ELECTRONS (NUM_MOLECULES * NUM_ELECTRONS_PER_MOLECULE)
#define TRAIL_LEN 120
#define NUCLEUS_RADIUS 0.6f
#define MOLECULE_SEPARATION 8.0f

// Physics constants
#define ATTRACTION_STRENGTH 0.8f
#define REPULSION_STRENGTH 1.2f
#define INTERACTION_RANGE 12.0f
#define MOLECULAR_SPEED 0.3f

// Visual quality settings for perfect spheres
#define SPHERE_QUALITY 36  
#define GLOW_INTENSITY 0.35f
#define TRAIL_ALPHA_START 0.8f
#define TRAIL_ALPHA_END 0.0f
#define TRAIL_FADE_SPEED 0.012f

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
    int molecule_id; // Which molecule this electron belongs to
} Electron;

typedef struct {
    Vec3 position;
    Vec3 velocity;
    Vec3 force;
    float charge; // Positive for nucleus
    Vec3 color;
    float energy_level;
} Molecule;

// Global state
TrailPoint electronTrail[TOTAL_ELECTRONS][TRAIL_LEN];
int trailIndex[TOTAL_ELECTRONS] = {0};
Electron electrons[TOTAL_ELECTRONS];
Molecule molecules[NUM_MOLECULES];
float global_time = 0.0f;

// Initialize molecular system with interactions
void initializeMolecularSystem() {
    srand(time(NULL));
    
    // Define molecule colors
    Vec3 molecule_colors[NUM_MOLECULES] = {
        {1.0f, 0.4f, 0.2f}, // Orange
        {0.2f, 0.8f, 0.4f}, // Green  
        {0.4f, 0.3f, 1.0f}, // Blue
        {0.9f, 0.7f, 0.2f}  // Yellow
    };
    
    // Initialize molecules in different positions
    for(int i = 0; i < NUM_MOLECULES; i++) {
        float angle = (2.0f * M_PI * i) / NUM_MOLECULES;
        molecules[i].position.x = MOLECULE_SEPARATION * cosf(angle);
        molecules[i].position.y = (i % 2 == 0) ? 2.0f : -2.0f;
        molecules[i].position.z = MOLECULE_SEPARATION * sinf(angle);
        
        molecules[i].velocity = (Vec3){0, 0, 0};
        molecules[i].force = (Vec3){0, 0, 0};
        molecules[i].charge = 6.0f; // Positive charge
        molecules[i].color = molecule_colors[i];
        molecules[i].energy_level = 1.0f + ((float)rand()/RAND_MAX) * 0.5f;
    }
    
    // Initialize electrons for each molecule
    Vec3 electron_colors[NUM_ELECTRONS_PER_MOLECULE] = {
        {0.3f, 0.9f, 1.0f}, // Cyan
        {1.0f, 0.4f, 0.6f}, // Pink
        {0.5f, 1.0f, 0.3f}  // Light Green
    };
    
    for(int mol = 0; mol < NUM_MOLECULES; mol++) {
        for(int e = 0; e < NUM_ELECTRONS_PER_MOLECULE; e++) {
            int electron_id = mol * NUM_ELECTRONS_PER_MOLECULE + e;
            
            electrons[electron_id].radius = 2.0f + e * 0.8f;
            electrons[electron_id].speed = 0.6f + ((float)rand()/RAND_MAX) * 0.8f;
            electrons[electron_id].angle_x = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
            electrons[electron_id].angle_y = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
            electrons[electron_id].angle_z = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
            electrons[electron_id].molecule_id = mol;
            electrons[electron_id].color = electron_colors[e];
        }
    }
    
    // Initialize all trail points
    for(int i = 0; i < TOTAL_ELECTRONS; i++) {
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

// Calculate distance between two points
float calculateDistance(Vec3 a, Vec3 b) {
    float dx = a.x - b.x;
    float dy = a.y - b.y; 
    float dz = a.z - b.z;
    return sqrtf(dx*dx + dy*dy + dz*dz);
}

// Calculate molecular forces between molecules
void calculateMolecularForces() {
    // Reset forces
    for(int i = 0; i < NUM_MOLECULES; i++) {
        molecules[i].force = (Vec3){0, 0, 0};
    }
    
    // Calculate forces between each pair of molecules
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            Vec3 pos1 = molecules[i].position;
            Vec3 pos2 = molecules[j].position;
            
            float distance = calculateDistance(pos1, pos2);
            
            if(distance < INTERACTION_RANGE && distance > 0.1f) {
                // Calculate direction vector
                Vec3 direction = {
                    (pos2.x - pos1.x) / distance,
                    (pos2.y - pos1.y) / distance,
                    (pos2.z - pos1.z) / distance
                };
                
                float force_magnitude;
                
                // Van der Waals-like forces: attraction at medium distance, repulsion at close distance
                if(distance < 4.0f) {
                    // Strong repulsion at close distance
                    force_magnitude = -REPULSION_STRENGTH / (distance * distance);
                } else {
                    // Weak attraction at medium distance
                    force_magnitude = ATTRACTION_STRENGTH / (distance * distance * distance);
                }
                
                // Apply forces (Newton's 3rd law)
                molecules[i].force.x -= direction.x * force_magnitude;
                molecules[i].force.y -= direction.y * force_magnitude;
                molecules[i].force.z -= direction.z * force_magnitude;
                
                molecules[j].force.x += direction.x * force_magnitude;
                molecules[j].force.y += direction.y * force_magnitude;
                molecules[j].force.z += direction.z * force_magnitude;
            }
        }
    }
}

// Update molecular physics
void updateMolecularPhysics(float dt) {
    calculateMolecularForces();
    
    for(int i = 0; i < NUM_MOLECULES; i++) {
        // Update velocity based on forces (F = ma, assuming mass = 1)
        molecules[i].velocity.x += molecules[i].force.x * dt;
        molecules[i].velocity.y += molecules[i].force.y * dt;
        molecules[i].velocity.z += molecules[i].force.z * dt;
        
        // Apply damping to prevent runaway motion
        molecules[i].velocity.x *= 0.98f;
        molecules[i].velocity.y *= 0.98f;
        molecules[i].velocity.z *= 0.98f;
        
        // Update position
        molecules[i].position.x += molecules[i].velocity.x * dt * MOLECULAR_SPEED;
        molecules[i].position.y += molecules[i].velocity.y * dt * MOLECULAR_SPEED;
        molecules[i].position.z += molecules[i].velocity.z * dt * MOLECULAR_SPEED;
        
        // Boundary conditions (keep molecules in view)
        float boundary = 15.0f;
        if(fabs(molecules[i].position.x) > boundary) molecules[i].velocity.x *= -0.5f;
        if(fabs(molecules[i].position.y) > boundary) molecules[i].velocity.y *= -0.5f;
        if(fabs(molecules[i].position.z) > boundary) molecules[i].velocity.z *= -0.5f;
    }
}

// Draw interaction lines between close molecules
void drawMolecularInteractions() {
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(2.0f);
    
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            float distance = calculateDistance(molecules[i].position, molecules[j].position);
            
            if(distance < INTERACTION_RANGE) {
                // Color based on interaction type
                float intensity = 1.0f - (distance / INTERACTION_RANGE);
                
                if(distance < 4.0f) {
                    // Repulsion - red lines
                    glColor4f(1.0f, 0.2f, 0.2f, intensity * 0.6f);
                } else {
                    // Attraction - blue lines
                    glColor4f(0.2f, 0.6f, 1.0f, intensity * 0.4f);
                }
                
                glBegin(GL_LINES);
                glVertex3f(molecules[i].position.x, molecules[i].position.y, molecules[i].position.z);
                glVertex3f(molecules[j].position.x, molecules[j].position.y, molecules[j].position.z);
                glEnd();
            }
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// Draw individual molecule nucleus
void drawMoleculeNucleus(int molecule_id) {
    Vec3 pos = molecules[molecule_id].position;
    Vec3 color = molecules[molecule_id].color;
    
    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    
    // Pulsing effect based on energy
    float pulse = 1.0f + 0.05f * sinf(global_time * 2.0f + molecule_id);
    glScalef(pulse, pulse, pulse);
    
    // Draw glowing aura
    drawSoftGlow(NUCLEUS_RADIUS, color.x, color.y, color.z);
    
    // Draw nucleus sphere
    drawPerfectSphere(NUCLEUS_RADIUS, color.x, color.y, color.z, 1.0f);
    
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

// Calculate 3D orbital positions relative to molecule center
Vec3 calculateElectronPosition(int electron_id, float time) {
    Electron* e = &electrons[electron_id];
    Vec3 molecule_pos = molecules[e->molecule_id].position;
    Vec3 local_pos = {0};
    
    // Update angles for 3D movement
    e->angle_x += time * e->speed * 0.8f;
    e->angle_y += time * e->speed * 0.6f;
    e->angle_z += time * e->speed * 0.4f;
    
    // Create 3D orbital motion with different patterns for each electron
    int orbital_type = electron_id % NUM_ELECTRONS_PER_MOLECULE;
    
    if(orbital_type == 0) {
        // Inner orbital - circular
        local_pos.x = e->radius * cosf(e->angle_x);
        local_pos.y = e->radius * sinf(e->angle_y) * 0.6f;
        local_pos.z = e->radius * sinf(e->angle_x);
    } else if(orbital_type == 1) {
        // Middle orbital - figure-8 pattern
        local_pos.x = e->radius * cosf(e->angle_x) * cosf(e->angle_z * 0.5f);
        local_pos.y = e->radius * sinf(e->angle_y) * 0.8f;
        local_pos.z = e->radius * sinf(e->angle_x) * sinf(e->angle_z * 0.3f);
    } else {
        // Outer orbital - complex 3D pattern
        local_pos.x = e->radius * cosf(e->angle_x) * sinf(e->angle_z * 0.7f);
        local_pos.y = e->radius * sinf(e->angle_y) * cosf(e->angle_z * 0.4f);
        local_pos.z = e->radius * sinf(e->angle_x) * cosf(e->angle_z * 0.6f);
    }
    
    // Add quantum uncertainty
    float variation = 0.08f;
    local_pos.x += (((float)rand()/RAND_MAX) - 0.5f) * variation;
    local_pos.y += (((float)rand()/RAND_MAX) - 0.5f) * variation;
    local_pos.z += (((float)rand()/RAND_MAX) - 0.5f) * variation;
    
    // Transform to world coordinates
    Vec3 world_pos = {
        molecule_pos.x + local_pos.x,
        molecule_pos.y + local_pos.y,
        molecule_pos.z + local_pos.z
    };
    
    return world_pos;
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

    // Initialize molecular system
    initializeMolecularSystem();

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
    
    printf("=== Molecular Simulation Controls ===\n");
    printf("ESC - Exit\n");
    printf("Arrow Keys - Rotate camera\n");
    printf("Up/Down Arrow - Zoom in/out\n");
    printf("Space - Reset view\n");
    printf("Watch molecules interact with each other!\n");
    printf("Blue lines = Attraction, Red lines = Repulsion\n\n");

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

        float dt = 0.016f; // ~60 FPS timing
        global_time += dt;
        
        // Update molecular physics
        updateMolecularPhysics(dt);

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
        gluLookAt(cam_x, 5.0, cam_z, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0);
        
        // Dark space background
        glClearColor(0.01f, 0.01f, 0.03f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw interaction lines between molecules
        drawMolecularInteractions();

        // Render all molecules
        for(int mol = 0; mol < NUM_MOLECULES; mol++) {
            drawMoleculeNucleus(mol);
        }

        // Render all electrons with 3D motion and trails
        for(int i = 0; i < TOTAL_ELECTRONS; i++) {
            // Calculate new position relative to molecule
            Vec3 pos = calculateElectronPosition(i, dt);
            
            // Update trail system
            updateTrail(i, pos);
            
            // Draw fading trail
            drawElectronTrail(i);

            // Draw perfect sphere electron with soft glow
            glPushMatrix();
            glTranslatef(pos.x, pos.y, pos.z);
            Vec3 color = electrons[i].color;
            
            // Draw soft glow around electron
            drawSoftGlow(0.06f, color.x, color.y, color.z);
            
            // Draw perfect electron sphere
            drawPerfectSphere(0.08f, color.x, color.y, color.z, 0.95f);
            
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
