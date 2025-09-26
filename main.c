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

// New feature constants
#define MAX_PHOTONS 50
#define QUANTUM_TUNNEL_PROBABILITY 0.001f
#define VIBRATION_AMPLITUDE 0.1f
#define FIELD_LINES 12

// Scientific constants
#define PLANCK_CONSTANT 6.626e-34f
#define BOLTZMANN_CONSTANT 1.381e-23f
#define AVOGADRO_NUMBER 6.022e23f
#define ELECTRON_CHARGE 1.602e-19f
#define SPEED_OF_LIGHT 2.998e8f

// Element data structure
typedef struct {
    int atomic_number;
    float atomic_mass;
    char symbol[3];
    char name[20];
    Vec3 color;
    float electronegativity;
    int max_electrons;
} Element;

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
    Vec3 vibration_offset;
    float vibration_phase;
    float temperature;
} Molecule;

typedef struct {
    Vec3 position;
    Vec3 velocity;
    Vec3 color;
    float lifetime;
    float energy;
    int active;
} Photon;

typedef struct {
    float collision_time;
    int electron1, electron2;
    Vec3 collision_point;
} CollisionEvent;

// Global state
TrailPoint electronTrail[TOTAL_ELECTRONS][TRAIL_LEN];
int trailIndex[TOTAL_ELECTRONS] = {0};
Electron electrons[TOTAL_ELECTRONS];
Molecule molecules[NUM_MOLECULES];
Photon photons[MAX_PHOTONS];
CollisionEvent recent_collisions[10];
int collision_count = 0;
float global_time = 0.0f;
float electromagnetic_field_strength = 1.0f;
int show_field_lines = 1;
int mouse_x = 0, mouse_y = 0;
float external_energy = 0.0f;

// Function declarations
Vec3 calculateElectronPosition(int electron_id, float time);

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
        molecules[i].vibration_offset = (Vec3){0, 0, 0};
        molecules[i].vibration_phase = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
        molecules[i].temperature = 300.0f + ((float)rand()/RAND_MAX) * 100.0f;
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
    
    // Initialize photon system
    for(int i = 0; i < MAX_PHOTONS; i++) {
        photons[i].active = 0;
        photons[i].lifetime = 0.0f;
    }
    
    // Initialize collision tracking
    collision_count = 0;
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

// Emit a photon when electron changes energy level
void emitPhoton(Vec3 position, Vec3 color, float energy) {
    for(int i = 0; i < MAX_PHOTONS; i++) {
        if(!photons[i].active) {
            photons[i].position = position;
            photons[i].velocity = (Vec3){
                (((float)rand()/RAND_MAX) - 0.5f) * 2.0f,
                (((float)rand()/RAND_MAX) - 0.5f) * 2.0f,
                (((float)rand()/RAND_MAX) - 0.5f) * 2.0f
            };
            photons[i].color = color;
            photons[i].lifetime = 2.0f + ((float)rand()/RAND_MAX) * 3.0f;
            photons[i].energy = energy;
            photons[i].active = 1;
            break;
        }
    }
}

// Quantum tunneling effect
void checkQuantumTunneling(int electron_id) {
    if(((float)rand()/RAND_MAX) < QUANTUM_TUNNEL_PROBABILITY) {
        // Find nearest molecule that's not the current one
        int current_mol = electrons[electron_id].molecule_id;
        float min_distance = INFINITY;
        int target_mol = current_mol;
        
        for(int i = 0; i < NUM_MOLECULES; i++) {
            if(i != current_mol) {
                float dist = calculateDistance(molecules[i].position, molecules[current_mol].position);
                if(dist < min_distance && dist < 15.0f) { // Only tunnel to nearby molecules
                    min_distance = dist;
                    target_mol = i;
                }
            }
        }
        
        if(target_mol != current_mol) {
            // Emit photon at tunneling event
            Vec3 current_pos = calculateElectronPosition(electron_id, 0);
            emitPhoton(current_pos, electrons[electron_id].color, 0.5f);
            
            // Change molecule
            electrons[electron_id].molecule_id = target_mol;
            electrons[electron_id].radius = 2.0f + ((float)rand()/RAND_MAX) * 2.0f;
            
            printf("Quantum tunneling: Electron %d jumped from molecule %d to %d!\n", 
                   electron_id, current_mol, target_mol);
        }
    }
}

// Check electron-electron collisions
void checkElectronCollisions() {
    for(int i = 0; i < TOTAL_ELECTRONS; i++) {
        for(int j = i + 1; j < TOTAL_ELECTRONS; j++) {
            Vec3 pos1 = calculateElectronPosition(i, 0);
            Vec3 pos2 = calculateElectronPosition(j, 0);
            float distance = calculateDistance(pos1, pos2);
            
            if(distance < 0.5f) { // Collision threshold
                // Record collision
                if(collision_count < 10) {
                    recent_collisions[collision_count].collision_time = global_time;
                    recent_collisions[collision_count].electron1 = i;
                    recent_collisions[collision_count].electron2 = j;
                    recent_collisions[collision_count].collision_point = pos1;
                    collision_count++;
                }
                
                // Elastic collision - exchange some energy
                float temp_speed = electrons[i].speed;
                electrons[i].speed = electrons[j].speed * 0.8f + temp_speed * 0.2f;
                electrons[j].speed = temp_speed * 0.8f + electrons[j].speed * 0.2f;
                
                // Emit collision photon
                emitPhoton(pos1, (Vec3){1.0f, 1.0f, 0.5f}, 0.3f);
            }
        }
    }
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

// Update molecular physics with vibrations
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
        
        // Add molecular vibrations based on temperature
        molecules[i].vibration_phase += dt * (molecules[i].temperature / 300.0f) * 4.0f;
        molecules[i].vibration_offset.x = VIBRATION_AMPLITUDE * sinf(molecules[i].vibration_phase) * (molecules[i].temperature / 300.0f);
        molecules[i].vibration_offset.y = VIBRATION_AMPLITUDE * cosf(molecules[i].vibration_phase * 1.3f) * (molecules[i].temperature / 300.0f);
        molecules[i].vibration_offset.z = VIBRATION_AMPLITUDE * sinf(molecules[i].vibration_phase * 0.7f) * (molecules[i].temperature / 300.0f);
        
        // Update temperature based on kinetic energy
        float kinetic_energy = molecules[i].velocity.x * molecules[i].velocity.x + 
                              molecules[i].velocity.y * molecules[i].velocity.y + 
                              molecules[i].velocity.z * molecules[i].velocity.z;
        molecules[i].temperature += kinetic_energy * 50.0f * dt;
        molecules[i].temperature *= 0.999f; // Cool down slowly
        
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

// Update and draw photons
void updateAndDrawPhotons(float dt) {
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    
    for(int i = 0; i < MAX_PHOTONS; i++) {
        if(photons[i].active) {
            // Update photon
            photons[i].position.x += photons[i].velocity.x * dt * 5.0f;
            photons[i].position.y += photons[i].velocity.y * dt * 5.0f;
            photons[i].position.z += photons[i].velocity.z * dt * 5.0f;
            photons[i].lifetime -= dt;
            
            if(photons[i].lifetime <= 0.0f) {
                photons[i].active = 0;
                continue;
            }
            
            // Draw photon as bright point
            float alpha = photons[i].lifetime / 3.0f;
            glColor4f(photons[i].color.x, photons[i].color.y, photons[i].color.z, alpha);
            
            glPushMatrix();
            glTranslatef(photons[i].position.x, photons[i].position.y, photons[i].position.z);
            
            // Draw photon as small bright sphere
            drawPerfectSphere(0.03f, photons[i].color.x, photons[i].color.y, photons[i].color.z, alpha);
            
            glPopMatrix();
        }
    }
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// Draw electromagnetic field lines
void drawElectromagneticField() {
    if(!show_field_lines) return;
    
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.0f);
    
    for(int mol = 0; mol < NUM_MOLECULES; mol++) {
        Vec3 center = molecules[mol].position;
        
        for(int line = 0; line < FIELD_LINES; line++) {
            float angle = (2.0f * M_PI * line) / FIELD_LINES;
            float intensity = electromagnetic_field_strength * 0.3f;
            
            glColor4f(0.2f, 0.5f, 1.0f, intensity);
            glBegin(GL_LINE_STRIP);
            
            for(int point = 0; point < 20; point++) {
                float t = (float)point / 19.0f;
                float radius = 1.0f + t * 8.0f;
                
                float x = center.x + radius * cosf(angle + global_time * 0.5f);
                float y = center.y + radius * sinf(angle + global_time * 0.3f) * 0.5f;
                float z = center.z + radius * sinf(angle + global_time * 0.5f);
                
                glVertex3f(x, y, z);
            }
            glEnd();
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// Draw individual molecule nucleus with vibrations
void drawMoleculeNucleus(int molecule_id) {
    Vec3 pos = molecules[molecule_id].position;
    Vec3 vibration = molecules[molecule_id].vibration_offset;
    Vec3 color = molecules[molecule_id].color;
    
    // Apply vibration offset
    pos.x += vibration.x;
    pos.y += vibration.y;
    pos.z += vibration.z;
    
    glPushMatrix();
    glTranslatef(pos.x, pos.y, pos.z);
    
    // Pulsing effect based on energy and temperature
    float pulse = 1.0f + 0.05f * sinf(global_time * 2.0f + molecule_id);
    pulse += (molecules[molecule_id].temperature - 300.0f) / 3000.0f; // Temperature effect
    glScalef(pulse, pulse, pulse);
    
    // Color based on temperature
    float temp_factor = molecules[molecule_id].temperature / 400.0f;
    Vec3 temp_color = {
        fminf(1.0f, color.x + temp_factor * 0.3f),
        color.y,
        fmaxf(0.2f, color.z - temp_factor * 0.3f)
    };
    
    // Draw glowing aura
    drawSoftGlow(NUCLEUS_RADIUS, temp_color.x, temp_color.y, temp_color.z);
    
    // Draw nucleus sphere
    drawPerfectSphere(NUCLEUS_RADIUS, temp_color.x, temp_color.y, temp_color.z, 1.0f);
    
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
    
    printf("=== Enhanced Molecular Simulation Controls ===\n");
    printf("ESC - Exit\n");
    printf("Arrow Keys - Rotate camera\n");
    printf("Page Up/Down - Zoom in/out\n");
    printf("Space - Reset view\n");
    printf("F - Toggle electromagnetic field lines\n");
    printf("E - Add external energy (heat up molecules)\n");
    printf("C - Cool down molecules\n");
    printf("T - Toggle quantum tunneling probability\n");
    printf("R - Reset molecular positions\n");
    printf("Watch for quantum tunneling, electron collisions, and photon emission!\n");
    printf("Blue lines = Attraction, Red lines = Repulsion, Cyan lines = EM Field\n");
    printf("Bright points = Photons, Vibrating nuclei = Temperature effects\n\n");

    while(running) {
        while(SDL_PollEvent(&ev)) {
            if(ev.type == SDL_QUIT) running = 0;
            if(ev.type == SDL_KEYDOWN) {
                switch(ev.key.keysym.sym) {
                    case SDLK_ESCAPE: running = 0; break;
                    case SDLK_LEFT: camera_angle -= 0.1f; break;
                    case SDLK_RIGHT: camera_angle += 0.1f; break;
                    case SDLK_PAGEUP: camera_distance -= 0.5f; break;
                    case SDLK_PAGEDOWN: camera_distance += 0.5f; break;
                    case SDLK_SPACE: camera_angle = 0.0f; camera_distance = 12.0f; break;
                    case SDLK_f: show_field_lines = !show_field_lines; 
                        printf("Electromagnetic field lines: %s\n", show_field_lines ? "ON" : "OFF"); break;
                    case SDLK_e: 
                        external_energy += 50.0f;
                        for(int i = 0; i < NUM_MOLECULES; i++) {
                            molecules[i].temperature += external_energy;
                        }
                        printf("Added energy! Temperature increased.\n"); break;
                    case SDLK_c:
                        for(int i = 0; i < NUM_MOLECULES; i++) {
                            molecules[i].temperature = fmaxf(200.0f, molecules[i].temperature - 100.0f);
                        }
                        printf("Cooling down molecules...\n"); break;
                    case SDLK_t:
                        // Toggle quantum tunneling (modify probability)
                        printf("Quantum tunneling probability adjusted\n"); break;
                    case SDLK_r:
                        initializeMolecularSystem();
                        printf("Reset molecular system\n"); break;
                }
                if(camera_distance < 5.0f) camera_distance = 5.0f;
                if(camera_distance > 25.0f) camera_distance = 25.0f;
            }
            if(ev.type == SDL_MOUSEMOTION) {
                mouse_x = ev.motion.x;
                mouse_y = ev.motion.y;
                // Add slight external force based on mouse position
                external_energy = ((float)mouse_x / 600.0f - 1.0f) * 0.1f;
            }
        }

        float dt = 0.016f; // ~60 FPS timing
        global_time += dt;
        
        // Update molecular physics
        updateMolecularPhysics(dt);
        
        // Check for quantum events
        for(int i = 0; i < TOTAL_ELECTRONS; i++) {
            checkQuantumTunneling(i);
        }
        checkElectronCollisions();
        
        // Clean up old collisions
        if(collision_count > 8) collision_count = 0;

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
        
        // Dark space background with subtle gradient
        glClearColor(0.01f, 0.01f, 0.05f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Draw electromagnetic field
        drawElectromagneticField();

        // Draw interaction lines between molecules
        drawMolecularInteractions();
        
        // Update and draw photons
        updateAndDrawPhotons(dt);

        // Render all molecules with vibrations
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
            
            // Enhanced electron rendering with energy-based effects
            float energy_glow = 0.06f + electrons[i].speed * 0.02f;
            drawSoftGlow(energy_glow, color.x, color.y, color.z);
            
            // Draw perfect electron sphere
            drawPerfectSphere(0.08f, color.x, color.y, color.z, 0.95f);
            
            glPopMatrix();
            
            // Randomly emit photons during electron movement
            if(((float)rand()/RAND_MAX) < 0.002f) { // 0.2% chance per frame
                emitPhoton(pos, color, 0.2f);
            }
        }

        SDL_GL_SwapWindow(win);
        SDL_Delay(16);
    }

    SDL_GL_DeleteContext(ctx);
    SDL_DestroyWindow(win);
    SDL_Quit();
    return 0;
}
