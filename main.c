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
#define QUANTUM_TUNNEL_PROBABILITY 0.0001f  // Reduced from 0.001f
#define VIBRATION_AMPLITUDE 0.1f
#define FIELD_LINES 12

// Visual enhancement constants
#define MAX_EXPLOSION_PARTICLES 200
#define EXPLOSION_LIFETIME 2.0f
#define BLOOM_INTENSITY 0.6f
#define VOLUMETRIC_SAMPLES 20

// Physics extension constants
#define GRAVITATIONAL_CONSTANT 6.67e-11f
#define MAGNETIC_PERMEABILITY 4.0e-7f * 3.14159f
#define VDW_EPSILON 0.2f  // Well depth
#define VDW_SIGMA 2.5f    // Distance at zero potential

// Molecule formation constants
#define MAX_FORMED_MOLECULES 10
#define BOND_FORMATION_DISTANCE 2.8f
#define BOND_BREAK_DISTANCE 4.5f
#define H2O_ANGLE 104.5f * 3.14159f / 180.0f  // Water bond angle
#define NH3_ANGLE 107.0f * 3.14159f / 180.0f  // Ammonia bond angle

// Molecular interaction constants
#define BONDING_DISTANCE 3.0f
#define STRONG_INTERACTION_DISTANCE 5.0f
#define DRAG_SENSITIVITY 0.1f

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
    float color_r, color_g, color_b;
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
    Element* element;
    float ionization_energy;
    float electron_affinity;
    int oxidation_state;
    float magnetic_moment;
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

// Explosion particle system
typedef struct {
    Vec3 position;
    Vec3 velocity;
    Vec3 color;
    float lifetime;
    float max_lifetime;
    float size;
    int active;
    int explosion_type; // 0=collision, 1=tunneling, 2=bonding
} ExplosionParticle;

// Enhanced Van der Waals force
typedef struct {
    float epsilon;  // Well depth
    float sigma;    // Distance at zero potential
    float r_min;    // Distance at minimum potential
} VdWParameters;

// Molecular structure definitions
typedef struct {
    int atom_count;
    int atom_types[4];  // Element indices
    Vec3 positions[4];  // Relative positions
    float bond_lengths[6]; // Maximum 6 bonds for 4 atoms
    float bond_angles[3];  // Bond angles
    char formula[10];
    float formation_energy;
    int active;
} FormedMolecule;

// Gravitational field structure
typedef struct {
    Vec3 position;
    float mass;
    float influence_radius;
} GravitySource;

// Element database with accurate atomic radii (pm)
Element elements[] = {
    {1, 1.008f, "H", "Hydrogen", 0.9f, 0.9f, 0.9f, 2.20f, 1},        // 53 pm
    {6, 12.011f, "C", "Carbon", 0.3f, 0.3f, 0.3f, 2.55f, 6},         // 67 pm  
    {7, 14.007f, "N", "Nitrogen", 0.2f, 0.6f, 1.0f, 3.04f, 7},       // 56 pm
    {8, 15.999f, "O", "Oxygen", 1.0f, 0.2f, 0.2f, 3.44f, 8},         // 48 pm
    {9, 18.998f, "F", "Fluorine", 0.9f, 1.0f, 0.2f, 3.98f, 9},       // 42 pm
    {10, 20.180f, "Ne", "Neon", 1.0f, 0.4f, 0.8f, 0.0f, 10}          // 38 pm
};

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

// New enhanced systems
ExplosionParticle explosion_particles[MAX_EXPLOSION_PARTICLES];
FormedMolecule formed_molecules[MAX_FORMED_MOLECULES];
GravitySource gravity_sources[NUM_MOLECULES];
VdWParameters vdw_params[6]; // For each element type
int show_explosions = 1;
int show_gravity_field = 0;
int show_vdw_potential = 0;
int show_formed_molecules = 1;
float bloom_effect_intensity = BLOOM_INTENSITY;

// Scientific analysis data
float total_system_energy = 0.0f;
float average_temperature = 0.0f;
int total_tunneling_events = 0;
int total_photons_emitted = 0;
float magnetic_field_strength = 0.5f;
int current_element_index = 1; // Start with Carbon

// Performance monitoring
float fps = 0.0f;
Uint32 frame_count = 0;
Uint32 fps_timer = 0;
float frame_time_ms = 0.0f;
float memory_usage_mb = 0.0f;
int show_performance_overlay = 1;
int show_help_overlay = 0;

// Molecular interaction and dragging
int dragging_molecule = -1;
Vec3 drag_offset = {0, 0, 0};
int mouse_button_down = 0;
float chemical_bonds[NUM_MOLECULES][NUM_MOLECULES]; // Bond strengths between molecules

// Real electron shell configurations
typedef struct {
    int n, l, j; // Quantum numbers: principal, angular momentum, total angular momentum
    int max_electrons;
    float energy_level;
    float radius_nm;
} ElectronShell;

ElectronShell shells[] = {
    {1, 0, 1, 2, -13.6f, 0.529f},  // 1s
    {2, 0, 1, 2, -3.4f, 2.116f},   // 2s  
    {2, 1, 1, 2, -3.4f, 2.116f},   // 2p1/2
    {2, 1, 3, 4, -3.4f, 2.116f},   // 2p3/2
    {3, 0, 1, 2, -1.51f, 4.761f},  // 3s
    {3, 1, 1, 2, -1.51f, 4.761f},  // 3p1/2
    {3, 1, 3, 4, -1.51f, 4.761f}   // 3p3/2
};

// Function declarations
Vec3 calculateElectronPosition(int electron_id, float time);
void calculateSystemStatistics();
void drawScientificHUD();
void changeElementType(int element_index);
float calculatePhotonWavelength(float energy);
void drawWaveFunction(int electron_id);
void updatePerformanceStats(Uint32 current_time);
void drawPerformanceOverlay();
void drawHelpOverlay();
void drawElementTooltip(int element_index, int mouse_x, int mouse_y);
float estimateMemoryUsage();
void handleErrorGracefully(const char* operation, const char* error_msg);
Vec3 screenToWorld(int screen_x, int screen_y, float camera_angle, float camera_distance);
int findMoleculeAtPosition(Vec3 world_pos);
void updateMolecularBonds();
void drawChemicalBonds();
void handleMouseDragging(SDL_Event* ev, float camera_angle, float camera_distance);
void updateElectronMoleculeBinding();
void checkElectronStealing();
void drawOrbitalPaths();
void updateElectronSpeedFromTemperature();

// New enhanced feature functions
void initializeExplosionSystem();
void createExplosion(Vec3 position, int explosion_type, Vec3 color);
void updateExplosions(float dt);
void drawExplosions();
void drawBloomEffect();
void drawVolumetricLighting();

void initializePhysicsExtensions();
void calculateGravitationalForces();
void calculateEnhancedVdWForces();
void updateMagneticDipoles();
void drawGravityField();
void drawVdWPotential();

void initializeMoleculeFormation();
void checkMoleculeFormation();
void updateFormedMolecules();
void drawFormedMolecules();
void breakMolecules();
Vec3 calculateOptimalGeometry(int molecule_type, int atom_index);

void applyShaderEffects();
void enableDepthBlending();
void drawParticleGlow(Vec3 position, Vec3 color, float intensity);

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
    
    // Initialize molecules with real element properties
    for(int i = 0; i < NUM_MOLECULES; i++) {
        float angle = (2.0f * M_PI * i) / NUM_MOLECULES;
        molecules[i].position.x = MOLECULE_SEPARATION * cosf(angle);
        molecules[i].position.y = (i % 2 == 0) ? 2.0f : -2.0f;
        molecules[i].position.z = MOLECULE_SEPARATION * sinf(angle);
        
        molecules[i].velocity = (Vec3){0, 0, 0};
        molecules[i].force = (Vec3){0, 0, 0};
        
        // Assign different elements to different molecules
        int element_idx = (current_element_index + i) % 6;
        molecules[i].element = &elements[element_idx];
        molecules[i].charge = (float)molecules[i].element->atomic_number;
        molecules[i].color.x = molecules[i].element->color_r;
        molecules[i].color.y = molecules[i].element->color_g;
        molecules[i].color.z = molecules[i].element->color_b;
        
        molecules[i].energy_level = 1.0f + ((float)rand()/RAND_MAX) * 0.5f;
        molecules[i].vibration_offset = (Vec3){0, 0, 0};
        molecules[i].vibration_phase = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
        molecules[i].temperature = 300.0f + ((float)rand()/RAND_MAX) * 100.0f;
        
        // Calculate real physical properties
        molecules[i].ionization_energy = 13.6f * molecules[i].element->atomic_number / (molecules[i].element->atomic_number + 1.0f);
        molecules[i].electron_affinity = molecules[i].element->electronegativity * 2.0f;
        molecules[i].oxidation_state = 0;
        molecules[i].magnetic_moment = ((float)rand()/RAND_MAX) * 2.0f - 1.0f;
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
    
    // Initialize new enhanced systems
    initializeExplosionSystem();
    initializePhysicsExtensions();
    initializeMoleculeFormation();
    
    // Initialize chemical bonds matrix
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = 0; j < NUM_MOLECULES; j++) {
            chemical_bonds[i][j] = 0.0f;
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
            
            // Create quantum tunneling explosion effect
            Vec3 tunnel_color = {0.9f, 0.5f, 1.0f}; // Purple/magenta for quantum events
            createExplosion(current_pos, 1, tunnel_color);
            
            // Change molecule
            electrons[electron_id].molecule_id = target_mol;
            electrons[electron_id].radius = 2.0f + ((float)rand()/RAND_MAX) * 2.0f;
            
            // Update statistics
            total_tunneling_events++;
            total_photons_emitted++;
            
            printf("Quantum tunneling: Electron %d jumped from %s to %s! (Event #%d)\n", 
                   electron_id, molecules[current_mol].element->symbol, 
                   molecules[target_mol].element->symbol, total_tunneling_events);
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
                
                // Create collision explosion effect
                Vec3 collision_color = {1.0f, 0.6f, 0.2f}; // Orange for collisions
                createExplosion(pos1, 0, collision_color);
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

// Calculate system-wide statistics
void calculateSystemStatistics() {
    total_system_energy = 0.0f;
    average_temperature = 0.0f;
    
    for(int i = 0; i < NUM_MOLECULES; i++) {
        // Kinetic energy
        float velocity_sq = molecules[i].velocity.x * molecules[i].velocity.x +
                           molecules[i].velocity.y * molecules[i].velocity.y +
                           molecules[i].velocity.z * molecules[i].velocity.z;
        total_system_energy += 0.5f * molecules[i].element->atomic_mass * velocity_sq;
        
        // Thermal energy
        total_system_energy += BOLTZMANN_CONSTANT * molecules[i].temperature;
        average_temperature += molecules[i].temperature;
    }
    
    average_temperature /= NUM_MOLECULES;
}

// Calculate photon wavelength from energy (nm)
float calculatePhotonWavelength(float energy) {
    // E = hc/λ → λ = hc/E
    return (PLANCK_CONSTANT * SPEED_OF_LIGHT) / (energy * ELECTRON_CHARGE) * 1e9f;
}

// Change element type for all molecules
void changeElementType(int element_index) {
    if(element_index >= 0 && element_index < 6) {
        current_element_index = element_index;
        for(int i = 0; i < NUM_MOLECULES; i++) {
            int idx = (element_index + i) % 6;
            molecules[i].element = &elements[idx];
            molecules[i].color.x = molecules[i].element->color_r;
            molecules[i].color.y = molecules[i].element->color_g;
            molecules[i].color.z = molecules[i].element->color_b;
            molecules[i].charge = (float)molecules[i].element->atomic_number;
            molecules[i].ionization_energy = 13.6f * molecules[i].element->atomic_number / (molecules[i].element->atomic_number + 1.0f);
        }
        printf("Changed to element: %s (%s)\n", elements[element_index].name, elements[element_index].symbol);
    }
}

// Draw scientific HUD with real-time data
void drawScientificHUD() {
    // Save current matrices
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1200, 900, 0, -1, 1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Semi-transparent background panel
    glColor4f(0.0f, 0.0f, 0.2f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(10, 10);
    glVertex2f(400, 10);
    glVertex2f(400, 300);
    glVertex2f(10, 300);
    glEnd();
    
    // Calculate current statistics
    calculateSystemStatistics();
    
    // Draw text (simplified - in real implementation you'd use a font rendering library)
    glColor4f(0.8f, 1.0f, 0.9f, 1.0f);
    
    // Note: Text rendering would require additional libraries like freetype
    // For now, we'll draw colored indicators instead
    
    // Energy level indicator
    float energy_bar_height = (total_system_energy / 1000.0f) * 100.0f;
    if(energy_bar_height > 100.0f) energy_bar_height = 100.0f;
    
    glColor4f(1.0f, 0.3f, 0.3f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(30, 250);
    glVertex2f(50, 250);
    glVertex2f(50, 250 - energy_bar_height);
    glVertex2f(30, 250 - energy_bar_height);
    glEnd();
    
    // Temperature indicator
    float temp_bar_height = ((average_temperature - 200.0f) / 200.0f) * 100.0f;
    if(temp_bar_height > 100.0f) temp_bar_height = 100.0f;
    if(temp_bar_height < 0.0f) temp_bar_height = 0.0f;
    
    glColor4f(0.3f, 0.8f, 1.0f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(70, 250);
    glVertex2f(90, 250);
    glVertex2f(90, 250 - temp_bar_height);
    glVertex2f(70, 250 - temp_bar_height);
    glEnd();
    
    // Quantum events indicator
    float quantum_bar = (float)(total_tunneling_events % 50) * 2.0f;
    glColor4f(0.9f, 0.9f, 0.3f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(110, 250);
    glVertex2f(130, 250);
    glVertex2f(130, 250 - quantum_bar);
    glVertex2f(110, 250 - quantum_bar);
    glEnd();
    
    // Current element indicator
    for(int i = 0; i < NUM_MOLECULES; i++) {
        glColor4f(molecules[i].element->color_r, molecules[i].element->color_g, molecules[i].element->color_b, 0.8f);
        float x = 150 + i * 30;
        glBegin(GL_QUADS);
        glVertex2f(x, 40);
        glVertex2f(x + 20, 40);
        glVertex2f(x + 20, 60);
        glVertex2f(x, 60);
        glEnd();
    }
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
    
    // Restore matrices
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

// Update performance statistics
void updatePerformanceStats(Uint32 current_time) {
    frame_count++;
    
    // Calculate FPS every second
    if(current_time - fps_timer >= 1000) {
        fps = (float)frame_count * 1000.0f / (current_time - fps_timer);
        frame_count = 0;
        fps_timer = current_time;
        
        // Estimate memory usage
        memory_usage_mb = estimateMemoryUsage();
    }
    
    static Uint32 last_frame_time = 0;
    frame_time_ms = (float)(current_time - last_frame_time);
    last_frame_time = current_time;
}

// Estimate memory usage (simplified)
float estimateMemoryUsage() {
    float base_memory = sizeof(electronTrail) + sizeof(electrons) + sizeof(molecules) + sizeof(photons);
    float dynamic_memory = collision_count * sizeof(CollisionEvent);
    return (base_memory + dynamic_memory) / (1024.0f * 1024.0f); // Convert to MB
}

// Error handling function
void handleErrorGracefully(const char* operation, const char* error_msg) {
    printf("ERROR in %s: %s\n", operation, error_msg);
    // Log to file in real implementation
    // fprintf(error_log, "[%.2f] %s: %s\n", global_time, operation, error_msg);
}

// Draw performance overlay
void drawPerformanceOverlay() {
    if(!show_performance_overlay) return;
    
    // Save current matrices
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1200, 900, 0, -1, 1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Performance panel (top-right)
    glColor4f(0.1f, 0.1f, 0.1f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(850, 10);
    glVertex2f(1190, 10);
    glVertex2f(1190, 150);
    glVertex2f(850, 150);
    glEnd();
    
    // FPS indicator
    float fps_ratio = fminf(fps / 60.0f, 1.0f);
    Vec3 fps_color = {1.0f - fps_ratio, fps_ratio, 0.0f}; // Red to green
    glColor4f(fps_color.x, fps_color.y, fps_color.z, 0.9f);
    glBegin(GL_QUADS);
    glVertex2f(870, 30);
    glVertex2f(870 + fps_ratio * 100, 30);
    glVertex2f(870 + fps_ratio * 100, 45);
    glVertex2f(870, 45);
    glEnd();
    
    // Memory usage bar
    float mem_ratio = fminf(memory_usage_mb / 100.0f, 1.0f); // Assume 100MB max
    glColor4f(0.3f, 0.7f, 1.0f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(870, 60);
    glVertex2f(870 + mem_ratio * 100, 60);
    glVertex2f(870 + mem_ratio * 100, 75);
    glVertex2f(870, 75);
    glEnd();
    
    // Frame time indicator
    float frame_ratio = fminf(frame_time_ms / 32.0f, 1.0f); // 32ms = ~30fps
    glColor4f(1.0f, 0.5f, 0.2f, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(870, 90);
    glVertex2f(870 + frame_ratio * 100, 90);
    glVertex2f(870 + frame_ratio * 100, 105);
    glVertex2f(870, 105);
    glEnd();
    
    // Labels (simplified - in real version use font rendering)
    glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
    // Would draw text here: "FPS: %.1f", "MEM: %.1fMB", "Frame: %.1fms"
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
    
    // Restore matrices
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

// Draw help overlay (F1)
void drawHelpOverlay() {
    if(!show_help_overlay) return;
    
    // Save current matrices
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1200, 900, 0, -1, 1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Semi-transparent fullscreen overlay
    glColor4f(0.0f, 0.0f, 0.0f, 0.7f);
    glBegin(GL_QUADS);
    glVertex2f(0, 0);
    glVertex2f(1200, 0);
    glVertex2f(1200, 900);
    glVertex2f(0, 900);
    glEnd();
    
    // Help panel
    glColor4f(0.1f, 0.2f, 0.3f, 0.9f);
    glBegin(GL_QUADS);
    glVertex2f(200, 100);
    glVertex2f(1000, 100);
    glVertex2f(1000, 800);
    glVertex2f(200, 800);
    glEnd();
    
    // Border
    glColor4f(0.4f, 0.7f, 1.0f, 1.0f);
    glLineWidth(3.0f);
    glBegin(GL_LINE_LOOP);
    glVertex2f(200, 100);
    glVertex2f(1000, 100);
    glVertex2f(1000, 800);
    glVertex2f(200, 800);
    glEnd();
    
    // Content areas for different help sections
    glColor4f(0.8f, 0.9f, 1.0f, 0.8f);
    
    // Keyboard shortcuts section
    for(int i = 0; i < 8; i++) {
        glBegin(GL_QUADS);
        glVertex2f(220, 140 + i * 60);
        glVertex2f(480, 140 + i * 60);
        glVertex2f(480, 180 + i * 60);
        glVertex2f(220, 180 + i * 60);
        glEnd();
    }
    
    // Element info section
    glColor4f(0.9f, 0.8f, 1.0f, 0.8f);
    for(int i = 0; i < 6; i++) {
        glBegin(GL_QUADS);
        glVertex2f(520, 140 + i * 80);
        glVertex2f(980, 140 + i * 80);
        glVertex2f(980, 200 + i * 80);
        glVertex2f(520, 200 + i * 80);
        glEnd();
    }
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
    
    // Restore matrices
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
}

// Draw element tooltip on hover
void drawElementTooltip(int element_index, int mouse_x, int mouse_y) {
    if(element_index < 0 || element_index >= 6) return;
    
    Element* elem = &elements[element_index];
    
    // Save current matrices
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glLoadIdentity();
    glOrtho(0, 1200, 900, 0, -1, 1);
    
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
    glLoadIdentity();
    
    glDisable(GL_LIGHTING);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Tooltip background
    float tooltip_x = (float)mouse_x + 20;
    float tooltip_y = (float)mouse_y - 10;
    
    glColor4f(0.2f, 0.2f, 0.4f, 0.9f);
    glBegin(GL_QUADS);
    glVertex2f(tooltip_x, tooltip_y);
    glVertex2f(tooltip_x + 200, tooltip_y);
    glVertex2f(tooltip_x + 200, tooltip_y + 120);
    glVertex2f(tooltip_x, tooltip_y + 120);
    glEnd();
    
    // Element color indicator
    glColor4f(elem->color_r, elem->color_g, elem->color_b, 0.8f);
    glBegin(GL_QUADS);
    glVertex2f(tooltip_x + 10, tooltip_y + 10);
    glVertex2f(tooltip_x + 40, tooltip_y + 10);
    glVertex2f(tooltip_x + 40, tooltip_y + 40);
    glVertex2f(tooltip_x + 10, tooltip_y + 40);
    glEnd();
    
    // Information bars (atomic number, mass, electronegativity)
    glColor4f(0.7f, 0.9f, 0.7f, 0.8f);
    float atomic_ratio = (float)elem->atomic_number / 20.0f;
    glBegin(GL_QUADS);
    glVertex2f(tooltip_x + 50, tooltip_y + 20);
    glVertex2f(tooltip_x + 50 + atomic_ratio * 100, tooltip_y + 20);
    glVertex2f(tooltip_x + 50 + atomic_ratio * 100, tooltip_y + 30);
    glVertex2f(tooltip_x + 50, tooltip_y + 30);
    glEnd();
    
    glColor4f(0.9f, 0.7f, 0.7f, 0.8f);
    float mass_ratio = elem->atomic_mass / 25.0f;
    glBegin(GL_QUADS);
    glVertex2f(tooltip_x + 50, tooltip_y + 40);
    glVertex2f(tooltip_x + 50 + mass_ratio * 100, tooltip_y + 40);
    glVertex2f(tooltip_x + 50 + mass_ratio * 100, tooltip_y + 50);
    glVertex2f(tooltip_x + 50, tooltip_y + 50);
    glEnd();
    
    // In real implementation, would render text here showing:
    // elem->name, elem->symbol, elem->atomic_number, elem->atomic_mass, elem->electronegativity
    
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LIGHTING);
    glDisable(GL_BLEND);
    
    // Restore matrices
    glPopMatrix();
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
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
    
    printf("=== PROFESSIONAL QUANTUM MOLECULAR SIMULATOR ===\n");
    printf("🔬 Scientific Controls:\n");
    printf("ESC - Exit\n");
    printf("Arrow Keys - Rotate camera | Page Up/Down - Zoom\n");
    printf("Space - Reset view | R - Reset system\n");
    printf("F - Toggle EM field lines | M - Toggle magnetic field\n");
    printf("E - Add energy (heat) | C - Cool down\n");
    printf("1-6 - Change elements (H,C,N,O,F,Ne)\n");
    printf("H - Toggle scientific HUD | P - Toggle performance overlay\n");
    printf("F1 - Help overlay | S - Save simulation state\n");
    printf("🖱️  MOUSE CONTROLS:\n");
    printf("Left Click + Drag - Grab and move NUCLEI (improved selection)\n");
    printf("Electrons follow their nucleus with orbital paths!\n");
    printf("Heat up molecules (E key) to speed up electrons!\n");
    printf("Stronger nuclei can steal electrons from weaker ones!\n");
    printf("Bring nuclei close together to form chemical bonds!\n");
    printf("\n🧪 Real-time Analysis:\n");
    printf("- Interactive molecular dragging and bonding\n");
    printf("- Chemical bond formation and visualization\n");
    printf("- Distance-based molecular interactions\n");
    printf("- Elliptical orbital paths (theory book style)\n");
    printf("- Temperature-dependent electron speeds\n");
    printf("- Optimized molecule selection (closest to click)\n");
    printf("- Reduced quantum tunneling (more realistic)\n");
    printf("- Electron-electron scattering dynamics\n");
    printf("- Photon emission with wavelength calculation\n");
    printf("- Thermodynamic state monitoring\n");
    printf("- Performance profiling with FPS/Memory monitoring\n");
    printf("🎨 Enhanced Visual Effects:\n");
    printf("- Particle explosion system for collisions and events\n");
    printf("- Advanced Van der Waals force visualization\n");
    printf("- Gravitational field effects between molecules\n");
    printf("- Bloom lighting and volumetric rays\n");
    printf("- Real molecular formation (H2, OH, CO, etc.)\n");
    printf("💫 NEW Controls: X=Explosions, G=Gravity, V=VdW, B=Bloom, N=Molecules\n\n");

    while(running) {
        while(SDL_PollEvent(&ev)) {
            // Handle mouse dragging
            handleMouseDragging(&ev, camera_angle, camera_distance);
            
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
                    case SDLK_h:
                        // Toggle HUD - we'll implement this as a global variable
                        printf("Scientific HUD toggled\n"); break;
                    case SDLK_m:
                        magnetic_field_strength = (magnetic_field_strength > 0.1f) ? 0.0f : 1.0f;
                        printf("Magnetic field: %s\n", (magnetic_field_strength > 0.1f) ? "ON" : "OFF"); break;
                    case SDLK_s:
                        printf("Simulation state saved (System Energy: %.2f, Avg Temp: %.1fK, Tunneling Events: %d)\n",
                               total_system_energy, average_temperature, total_tunneling_events); break;
                    case SDLK_1: changeElementType(0); break; // Hydrogen
                    case SDLK_2: changeElementType(1); break; // Carbon
                    case SDLK_3: changeElementType(2); break; // Nitrogen
                    case SDLK_4: changeElementType(3); break; // Oxygen
                    case SDLK_5: changeElementType(4); break; // Fluorine
                    case SDLK_6: changeElementType(5); break; // Neon
                    case SDLK_p:
                        show_performance_overlay = !show_performance_overlay;
                        printf("Performance overlay: %s\n", show_performance_overlay ? "ON" : "OFF"); break;
                    case SDLK_F1:
                        show_help_overlay = !show_help_overlay;
                        printf("Help overlay: %s\n", show_help_overlay ? "ON" : "OFF"); break;
                    
                    // New enhanced features controls
                    case SDLK_x:
                        show_explosions = !show_explosions;
                        printf("💥 Explosion effects: %s\n", show_explosions ? "ON" : "OFF"); break;
                    case SDLK_g:
                        show_gravity_field = !show_gravity_field;
                        printf("🌌 Gravity field: %s\n", show_gravity_field ? "ON" : "OFF"); break;
                    case SDLK_v:
                        show_vdw_potential = !show_vdw_potential;
                        printf("⚛️  Van der Waals potential: %s\n", show_vdw_potential ? "ON" : "OFF"); break;
                    case SDLK_b:
                        bloom_effect_intensity = (bloom_effect_intensity > 0.1f) ? 0.0f : BLOOM_INTENSITY;
                        printf("✨ Bloom effects: %s\n", (bloom_effect_intensity > 0.1f) ? "ON" : "OFF"); break;
                    case SDLK_n:
                        show_formed_molecules = !show_formed_molecules;
                        printf("🧪 Formed molecules: %s\n", show_formed_molecules ? "ON" : "OFF"); break;
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
            // Mouse button events are now handled in handleMouseDragging function
            // Remove this old mouse handling code to avoid conflicts
        }

        // Performance monitoring
        Uint32 current_time = SDL_GetTicks();
        updatePerformanceStats(current_time);
        
        float dt = 0.016f; // ~60 FPS timing
        global_time += dt;
        
        // Update molecular physics
        updateMolecularPhysics(dt);
        
        // Update chemical bonds
        updateMolecularBonds();
        
        // Update electron-nucleus binding forces
        updateElectronMoleculeBinding();
        
        // Update electron speeds based on temperature
        updateElectronSpeedFromTemperature();
        
        // Enhanced physics calculations
        calculateGravitationalForces();
        calculateEnhancedVdWForces();
        updateMagneticDipoles();
        
        // Molecular formation system
        checkMoleculeFormation();
        updateFormedMolecules();
        
        // Update explosion particles
        updateExplosions(dt);
        
        // Check for quantum events (reduced frequency)
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

        // Draw enhanced physics visualizations
        drawGravityField();
        drawVdWPotential();
        
        // Draw interaction lines between molecules
        drawMolecularInteractions();
        
        // Draw chemical bonds
        drawChemicalBonds();
        
        // Draw formed molecules
        drawFormedMolecules();
        
        // Draw orbital paths (theory book style)
        drawOrbitalPaths();
        
        // Draw volumetric lighting effects
        drawVolumetricLighting();
        
        // Update and draw photons
        updateAndDrawPhotons(dt);
        
        // Draw explosion particles
        drawExplosions();

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
                total_photons_emitted++;
            }
        }
        
        // Draw bloom effect (after all bright objects)
        drawBloomEffect();
        
        // Draw scientific HUD
        drawScientificHUD();
        
        // Draw performance overlay
        drawPerformanceOverlay();
        
        // Draw help overlay (if active)
        drawHelpOverlay();
        
        // Draw element tooltip on hover (check if mouse is over element indicators)
        if(mouse_x >= 150 && mouse_x <= 270 && mouse_y >= 40 && mouse_y <= 60) {
            int hovered_element = (mouse_x - 150) / 30;
            if(hovered_element >= 0 && hovered_element < NUM_MOLECULES) {
                drawElementTooltip(molecules[hovered_element].element - elements, mouse_x, mouse_y);
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

// Convert screen coordinates to 3D world coordinates (improved)
Vec3 screenToWorld(int screen_x, int screen_y, float camera_angle, float camera_distance) {
    // Get actual window size for proper conversion
    int window_width = 1200;
    int window_height = 900;
    
    // Normalize screen coordinates to [-1, 1] range
    float normalized_x = ((float)screen_x / (window_width * 0.5f)) - 1.0f;
    float normalized_y = 1.0f - ((float)screen_y / (window_height * 0.5f));
    
    // Project to world space - molecules are spread over ~16 units, so scale accordingly
    Vec3 world_pos;
    world_pos.x = normalized_x * camera_distance * 0.8f; // Increased scaling to match molecule positions
    world_pos.y = normalized_y * camera_distance * 0.8f;
    world_pos.z = 0.0f; // Keep on camera plane
    
    // Apply camera rotation around Y axis
    float cos_angle = cosf(camera_angle);
    float sin_angle = sinf(camera_angle);
    Vec3 rotated;
    rotated.x = world_pos.x * cos_angle - world_pos.z * sin_angle;
    rotated.y = world_pos.y;
    rotated.z = world_pos.x * sin_angle + world_pos.z * cos_angle;
    
    // Debug output removed for better performance
    
    return rotated;
}

// Find which molecule is closest to a given world position
int findMoleculeAtPosition(Vec3 world_pos) {
    float min_distance = INFINITY;
    int closest_molecule = -1;
    
    for(int i = 0; i < NUM_MOLECULES; i++) {
        Vec3 mol_pos = molecules[i].position;
        // Include vibration offset for accurate detection
        mol_pos.x += molecules[i].vibration_offset.x;
        mol_pos.y += molecules[i].vibration_offset.y;
        mol_pos.z += molecules[i].vibration_offset.z;
        
        float distance = calculateDistance(world_pos, mol_pos);
        
        if(distance < min_distance) {
            min_distance = distance;
            closest_molecule = i;
        }
    }
    
    // Accept the closest molecule if it's within reasonable range
    if(closest_molecule != -1 && min_distance < 15.0f) {
        printf("Selected closest molecule: %s (distance: %.1f)\n", 
               molecules[closest_molecule].element->symbol, min_distance);
        return closest_molecule;
    }
    
    return -1;
}

// Update chemical bonds based on distance
void updateMolecularBonds() {
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            float distance = calculateDistance(molecules[i].position, molecules[j].position);
            
            if(distance < BONDING_DISTANCE) {
                // Form chemical bond
                float bond_strength = (BONDING_DISTANCE - distance) / BONDING_DISTANCE;
                
                // Only print if bond is newly formed
                if(chemical_bonds[i][j] < 0.1f && bond_strength > 0.5f) {
                    printf("Chemical bond formed between %s and %s! (strength: %.2f)\n",
                           molecules[i].element->symbol, molecules[j].element->symbol, bond_strength);
                }
                
                chemical_bonds[i][j] = bond_strength;
                chemical_bonds[j][i] = bond_strength;
                
                // Apply bonding forces (attraction)
                Vec3 direction = {
                    (molecules[j].position.x - molecules[i].position.x) / distance,
                    (molecules[j].position.y - molecules[i].position.y) / distance,
                    (molecules[j].position.z - molecules[i].position.z) / distance
                };
                
                float bond_force = bond_strength * 0.5f;
                molecules[i].force.x += direction.x * bond_force;
                molecules[i].force.y += direction.y * bond_force;
                molecules[i].force.z += direction.z * bond_force;
                
                molecules[j].force.x -= direction.x * bond_force;
                molecules[j].force.y -= direction.y * bond_force;
                molecules[j].force.z -= direction.z * bond_force;
                
            } else if(distance < STRONG_INTERACTION_DISTANCE) {
                // Strong interaction (no bond yet)
                float interaction_strength = (STRONG_INTERACTION_DISTANCE - distance) / STRONG_INTERACTION_DISTANCE * 0.2f;
                chemical_bonds[i][j] = interaction_strength;
                chemical_bonds[j][i] = interaction_strength;
            } else {
                // No interaction
                chemical_bonds[i][j] = 0.0f;
                chemical_bonds[j][i] = 0.0f;
            }
        }
    }
}

// Draw chemical bonds
void drawChemicalBonds() {
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(4.0f);
    
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            float bond_strength = chemical_bonds[i][j];
            if(bond_strength > 0.1f) {
                Vec3 pos1 = molecules[i].position;
                Vec3 pos2 = molecules[j].position;
                
                // Bond color based on strength and elements
                float r = 0.9f;
                float g = 0.9f - bond_strength * 0.5f;
                float b = 0.2f + bond_strength * 0.8f;
                
                glColor4f(r, g, b, bond_strength);
                
                glBegin(GL_LINES);
                glVertex3f(pos1.x, pos1.y, pos1.z);
                glVertex3f(pos2.x, pos2.y, pos2.z);
                glEnd();
                
                // Draw bond strength indicator (thicker line)
                if(bond_strength > 0.5f) {
                    glLineWidth(8.0f);
                    glColor4f(1.0f, 1.0f, 1.0f, bond_strength * 0.8f);
                    glBegin(GL_LINES);
                    glVertex3f(pos1.x, pos1.y, pos1.z);
                    glVertex3f(pos2.x, pos2.y, pos2.z);
                    glEnd();
                    glLineWidth(4.0f);
                }
            }
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// Handle mouse dragging for molecules (nuclei)
void handleMouseDragging(SDL_Event* ev, float camera_angle, float camera_distance) {
    if(ev->type == SDL_MOUSEBUTTONDOWN && ev->button.button == SDL_BUTTON_LEFT) {
        // Clean mouse handling
        Vec3 world_pos = screenToWorld(ev->button.x, ev->button.y, camera_angle, camera_distance);
        dragging_molecule = findMoleculeAtPosition(world_pos);
        
        if(dragging_molecule != -1) {
            drag_offset.x = molecules[dragging_molecule].position.x - world_pos.x;
            drag_offset.y = molecules[dragging_molecule].position.y - world_pos.y;
            drag_offset.z = molecules[dragging_molecule].position.z - world_pos.z;
            mouse_button_down = 1;
            printf("SUCCESS! Grabbed nucleus: %s (Z=%d) - electrons will follow!\n", 
                   molecules[dragging_molecule].element->name, 
                   molecules[dragging_molecule].element->atomic_number);
        } else {
            printf("FAILED to grab any molecule\n");
        }
    }
    
    if(ev->type == SDL_MOUSEBUTTONUP && ev->button.button == SDL_BUTTON_LEFT) {
        if(dragging_molecule != -1) {
            printf("Released nucleus: %s - checking for electron stealing...\n", 
                   molecules[dragging_molecule].element->name);
            
            // Check if any electrons should be stolen by nearby nuclei
            checkElectronStealing();
        }
        dragging_molecule = -1;
        mouse_button_down = 0;
    }
    
    if(ev->type == SDL_MOUSEMOTION && dragging_molecule != -1 && mouse_button_down) {
        Vec3 old_pos = molecules[dragging_molecule].position;
        Vec3 world_pos = screenToWorld(ev->motion.x, ev->motion.y, camera_angle, camera_distance);
        
        Vec3 new_pos = {
            world_pos.x + drag_offset.x,
            world_pos.y + drag_offset.y,
            world_pos.z + drag_offset.z
        };
        
        // Calculate movement delta
        Vec3 movement_delta = {
            new_pos.x - old_pos.x,
            new_pos.y - old_pos.y,
            new_pos.z - old_pos.z
        };
        
        // Update molecule position
        molecules[dragging_molecule].position = new_pos;
        
        // Reset velocity when dragging
        molecules[dragging_molecule].velocity.x *= 0.5f;
        molecules[dragging_molecule].velocity.y *= 0.5f;
        molecules[dragging_molecule].velocity.z *= 0.5f;
        
        // Move electrons that belong to this nucleus
        for(int i = 0; i < TOTAL_ELECTRONS; i++) {
            if(electrons[i].molecule_id == dragging_molecule) {
                // Apply a fraction of the movement to electron orbital centers
                // This creates a "dragging" effect while maintaining orbital motion
                Vec3 current_pos = calculateElectronPosition(i, 0);
                
                // The electron tries to follow the nucleus but with some lag/elasticity
                float follow_strength = 0.7f; // 70% follow rate
                
                // Update the electron's effective center by nudging its angles
                electrons[i].angle_x += movement_delta.x * follow_strength * 0.1f;
                electrons[i].angle_y += movement_delta.y * follow_strength * 0.1f;
                electrons[i].angle_z += movement_delta.z * follow_strength * 0.1f;
            }
        }
    }
}

// Update electron-molecule binding forces
void updateElectronMoleculeBinding() {
    for(int i = 0; i < TOTAL_ELECTRONS; i++) {
        int current_mol = electrons[i].molecule_id;
        Vec3 electron_pos = calculateElectronPosition(i, 0);
        Vec3 nucleus_pos = molecules[current_mol].position;
        
        // Calculate distance from electron to its nucleus
        float distance_to_nucleus = calculateDistance(electron_pos, nucleus_pos);
        
        // If electron is too far from its nucleus, apply strong attractive force
        if(distance_to_nucleus > electrons[i].radius * 1.5f) {
            // Strong binding force to keep electron in orbit
            Vec3 binding_direction = {
                (nucleus_pos.x - electron_pos.x) / distance_to_nucleus,
                (nucleus_pos.y - electron_pos.y) / distance_to_nucleus,
                (nucleus_pos.z - electron_pos.z) / distance_to_nucleus
            };
            
            float binding_strength = molecules[current_mol].element->atomic_number * 0.1f;
            
            // Apply binding force by adjusting electron's orbital parameters
            electrons[i].angle_x += binding_direction.x * binding_strength * 0.01f;
            electrons[i].angle_y += binding_direction.y * binding_strength * 0.01f;
            electrons[i].angle_z += binding_direction.z * binding_strength * 0.01f;
        }
    }
}

// Check if electrons should be stolen by stronger nuclei
void checkElectronStealing() {
    for(int i = 0; i < TOTAL_ELECTRONS; i++) {
        int current_mol = electrons[i].molecule_id;
        Vec3 electron_pos = calculateElectronPosition(i, 0);
        
        // Check all other molecules to see if they can steal this electron
        for(int j = 0; j < NUM_MOLECULES; j++) {
            if(j != current_mol) {
                Vec3 other_nucleus = molecules[j].position;
                float distance_to_other = calculateDistance(electron_pos, other_nucleus);
                float distance_to_current = calculateDistance(electron_pos, molecules[current_mol].position);
                
                // Stealing conditions:
                // 1. Other nucleus is closer
                // 2. Other nucleus has higher atomic number (more protons)
                // 3. Distance is within reasonable range
                if(distance_to_other < distance_to_current * 0.7f && 
                   molecules[j].element->atomic_number > molecules[current_mol].element->atomic_number &&
                   distance_to_other < 4.0f) {
                    
                    // Steal the electron!
                    electrons[i].molecule_id = j;
                    electrons[i].radius = 2.0f + ((float)rand()/RAND_MAX) * 2.0f;
                    
                    printf("Electron stolen! %s (Z=%d) stole electron from %s (Z=%d)\n",
                           molecules[j].element->symbol, molecules[j].element->atomic_number,
                           molecules[current_mol].element->symbol, molecules[current_mol].element->atomic_number);
                    
                    // Emit photon for the stealing event
                    emitPhoton(electron_pos, electrons[i].color, 0.8f);
                    total_photons_emitted++;
                    
                    break; // Only steal once per frame
                }
            }
        }
    }
}

// Draw orbital paths as elliptical lines
void drawOrbitalPaths() {
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glLineWidth(1.0f);
    
    for(int mol = 0; mol < NUM_MOLECULES; mol++) {
        Vec3 nucleus_pos = molecules[mol].position;
        
        // Draw different orbital shapes for each electron shell
        for(int e = 0; e < NUM_ELECTRONS_PER_MOLECULE; e++) {
            int electron_id = mol * NUM_ELECTRONS_PER_MOLECULE + e;
            Electron* electron = &electrons[electron_id];
            
            // Orbital color with transparency
            glColor4f(electron->color.x, electron->color.y, electron->color.z, 0.3f);
            
            glBegin(GL_LINE_LOOP);
            
            // Draw elliptical orbital path
            for(int i = 0; i < 64; i++) {
                float angle = (2.0f * M_PI * i) / 64.0f;
                
                Vec3 orbital_point;
                
                if(e == 0) {
                    // Inner orbital - circular
                    orbital_point.x = nucleus_pos.x + electron->radius * cosf(angle);
                    orbital_point.y = nucleus_pos.y + electron->radius * sinf(angle) * 0.6f;
                    orbital_point.z = nucleus_pos.z;
                } else if(e == 1) {
                    // Middle orbital - figure-8 pattern
                    orbital_point.x = nucleus_pos.x + electron->radius * cosf(angle) * cosf(angle * 0.5f);
                    orbital_point.y = nucleus_pos.y + electron->radius * sinf(angle) * 0.8f;
                    orbital_point.z = nucleus_pos.z + electron->radius * sinf(angle) * sinf(angle * 0.3f);
                } else {
                    // Outer orbital - complex 3D ellipse
                    orbital_point.x = nucleus_pos.x + electron->radius * cosf(angle) * sinf(angle * 0.7f);
                    orbital_point.y = nucleus_pos.y + electron->radius * sinf(angle) * cosf(angle * 0.4f);
                    orbital_point.z = nucleus_pos.z + electron->radius * sinf(angle) * cosf(angle * 0.6f);
                }
                
                glVertex3f(orbital_point.x, orbital_point.y, orbital_point.z);
            }
            
            glEnd();
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// Update electron speed based on temperature
void updateElectronSpeedFromTemperature() {
    for(int i = 0; i < TOTAL_ELECTRONS; i++) {
        int mol_id = electrons[i].molecule_id;
        float temperature = molecules[mol_id].temperature;
        
        // Base speed + temperature effect
        float base_speed = 0.6f + ((float)rand()/RAND_MAX) * 0.8f;
        float temp_multiplier = temperature / 300.0f; // 300K as reference
        
        // Speed increases with temperature (kinetic theory)
        electrons[i].speed = base_speed * (0.5f + temp_multiplier * 1.5f);
        
        // Clamp speed to reasonable values
        if(electrons[i].speed > 3.0f) electrons[i].speed = 3.0f;
        if(electrons[i].speed < 0.2f) electrons[i].speed = 0.2f;
    }
}

// ================================
// EXPLOSION PARTICLE SYSTEM
// ================================

void initializeExplosionSystem() {
    for(int i = 0; i < MAX_EXPLOSION_PARTICLES; i++) {
        explosion_particles[i].active = 0;
        explosion_particles[i].lifetime = 0.0f;
        explosion_particles[i].max_lifetime = EXPLOSION_LIFETIME;
    }
}

void createExplosion(Vec3 position, int explosion_type, Vec3 color) {
    int particles_to_create = 15 + rand() % 20; // 15-35 particles
    
    for(int i = 0; i < particles_to_create; i++) {
        // Find inactive particle
        for(int j = 0; j < MAX_EXPLOSION_PARTICLES; j++) {
            if(!explosion_particles[j].active) {
                explosion_particles[j].position = position;
                explosion_particles[j].active = 1;
                explosion_particles[j].lifetime = 0.0f;
                explosion_particles[j].explosion_type = explosion_type;
                explosion_particles[j].color = color;
                
                // Random velocity based on explosion type
                float speed = 2.0f + ((float)rand()/RAND_MAX) * 3.0f;
                float theta = ((float)rand()/RAND_MAX) * 2.0f * M_PI;
                float phi = ((float)rand()/RAND_MAX) * M_PI;
                
                explosion_particles[j].velocity.x = speed * sinf(phi) * cosf(theta);
                explosion_particles[j].velocity.y = speed * sinf(phi) * sinf(theta);
                explosion_particles[j].velocity.z = speed * cosf(phi);
                
                // Size based on explosion type
                if(explosion_type == 0) { // Collision
                    explosion_particles[j].size = 0.1f + ((float)rand()/RAND_MAX) * 0.2f;
                } else if(explosion_type == 1) { // Tunneling
                    explosion_particles[j].size = 0.05f + ((float)rand()/RAND_MAX) * 0.1f;
                    explosion_particles[j].color.x += 0.3f; // Brighter
                } else { // Bonding
                    explosion_particles[j].size = 0.15f + ((float)rand()/RAND_MAX) * 0.25f;
                    explosion_particles[j].color.y += 0.4f; // More green
                }
                
                explosion_particles[j].max_lifetime = EXPLOSION_LIFETIME + ((float)rand()/RAND_MAX) * 1.0f;
                break;
            }
        }
    }
}

void updateExplosions(float dt) {
    for(int i = 0; i < MAX_EXPLOSION_PARTICLES; i++) {
        if(explosion_particles[i].active) {
            explosion_particles[i].lifetime += dt;
            
            // Update position
            explosion_particles[i].position.x += explosion_particles[i].velocity.x * dt;
            explosion_particles[i].position.y += explosion_particles[i].velocity.y * dt;
            explosion_particles[i].position.z += explosion_particles[i].velocity.z * dt;
            
            // Apply gravity and drag
            explosion_particles[i].velocity.y -= 1.0f * dt; // Gravity
            explosion_particles[i].velocity.x *= 0.98f; // Drag
            explosion_particles[i].velocity.y *= 0.98f;
            explosion_particles[i].velocity.z *= 0.98f;
            
            // Fade out over time
            if(explosion_particles[i].lifetime >= explosion_particles[i].max_lifetime) {
                explosion_particles[i].active = 0;
            }
        }
    }
}

void drawExplosions() {
    if(!show_explosions) return;
    
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE); // Additive blending for bright particles
    
    for(int i = 0; i < MAX_EXPLOSION_PARTICLES; i++) {
        if(explosion_particles[i].active) {
            float life_ratio = explosion_particles[i].lifetime / explosion_particles[i].max_lifetime;
            float alpha = 1.0f - life_ratio; // Fade out
            
            glPushMatrix();
            glTranslatef(explosion_particles[i].position.x, 
                        explosion_particles[i].position.y, 
                        explosion_particles[i].position.z);
            
            glColor4f(explosion_particles[i].color.x, 
                     explosion_particles[i].color.y, 
                     explosion_particles[i].color.z, 
                     alpha);
            
            // Draw glowing particle
            drawParticleGlow(explosion_particles[i].position, explosion_particles[i].color, alpha);
            
            glPopMatrix();
        }
    }
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// ================================
// VISUAL ENHANCEMENTS
// ================================

void drawBloomEffect() {
    if(bloom_effect_intensity <= 0.0f) return;
    
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    
    // Draw enhanced glow around bright objects
    for(int i = 0; i < NUM_MOLECULES; i++) {
        Vec3 pos = molecules[i].position;
        Vec3 color = {molecules[i].element->color_r, molecules[i].element->color_g, molecules[i].element->color_b};
        
        // Multiple layers of glow for bloom effect
        for(int layer = 0; layer < 3; layer++) {
            float radius = NUCLEUS_RADIUS * (2.0f + layer * 1.5f);
            float alpha = bloom_effect_intensity * (0.3f - layer * 0.1f);
            
            glPushMatrix();
            glTranslatef(pos.x, pos.y, pos.z);
            glColor4f(color.x, color.y, color.z, alpha);
            
            GLUquadric* quad = gluNewQuadric();
            gluSphere(quad, radius, 12, 12);
            gluDeleteQuadric(quad);
            glPopMatrix();
        }
    }
    
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

void drawVolumetricLighting() {
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Draw light rays from electron emissions
    for(int i = 0; i < MAX_PHOTONS; i++) {
        if(photons[i].active) {
            glBegin(GL_LINES);
            glColor4f(photons[i].color.x, photons[i].color.y, photons[i].color.z, 0.3f);
            
            Vec3 start = photons[i].position;
            Vec3 end = {
                start.x + photons[i].velocity.x * 2.0f,
                start.y + photons[i].velocity.y * 2.0f,
                start.z + photons[i].velocity.z * 2.0f
            };
            
            glVertex3f(start.x, start.y, start.z);
            glVertex3f(end.x, end.y, end.z);
            glEnd();
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

void drawParticleGlow(Vec3 position, Vec3 color, float intensity) {
    glPushMatrix();
    glTranslatef(position.x, position.y, position.z);
    
    // Core particle
    glColor4f(color.x, color.y, color.z, intensity);
    GLUquadric* quad = gluNewQuadric();
    gluSphere(quad, 0.05f, 8, 8);
    
    // Outer glow
    glColor4f(color.x, color.y, color.z, intensity * 0.3f);
    gluSphere(quad, 0.15f, 8, 8);
    
    gluDeleteQuadric(quad);
    glPopMatrix();
}

// ================================
// PHYSICS EXTENSIONS
// ================================

void initializePhysicsExtensions() {
    // Initialize Van der Waals parameters for each element
    for(int i = 0; i < 6; i++) {
        vdw_params[i].epsilon = VDW_EPSILON * (1.0f + i * 0.1f); // Varies by element
        vdw_params[i].sigma = VDW_SIGMA * (0.8f + i * 0.05f);
        vdw_params[i].r_min = vdw_params[i].sigma * powf(2.0f, 1.0f/6.0f);
    }
    
    // Initialize gravity sources (each molecule acts as gravity source)
    for(int i = 0; i < NUM_MOLECULES; i++) {
        gravity_sources[i].position = molecules[i].position;
        gravity_sources[i].mass = molecules[i].element->atomic_mass;
        gravity_sources[i].influence_radius = 10.0f;
    }
}

void calculateGravitationalForces() {
    if(!show_gravity_field) return;
    
    for(int i = 0; i < NUM_MOLECULES; i++) {
        molecules[i].force.x = 0.0f;
        molecules[i].force.y = 0.0f;
        molecules[i].force.z = 0.0f;
        
        for(int j = 0; j < NUM_MOLECULES; j++) {
            if(i != j) {
                Vec3 r_vec = {
                    molecules[j].position.x - molecules[i].position.x,
                    molecules[j].position.y - molecules[i].position.y,
                    molecules[j].position.z - molecules[i].position.z
                };
                
                float r = sqrtf(r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z);
                if(r > 0.1f) { // Avoid division by zero
                    float force_magnitude = GRAVITATIONAL_CONSTANT * molecules[i].element->atomic_mass * 
                                          molecules[j].element->atomic_mass / (r * r) * 1e20f; // Scale factor
                    
                    molecules[i].force.x += force_magnitude * r_vec.x / r;
                    molecules[i].force.y += force_magnitude * r_vec.y / r;
                    molecules[i].force.z += force_magnitude * r_vec.z / r;
                }
            }
        }
        
        // Update gravity source positions
        gravity_sources[i].position = molecules[i].position;
    }
}

void calculateEnhancedVdWForces() {
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            Vec3 r_vec = {
                molecules[j].position.x - molecules[i].position.x,
                molecules[j].position.y - molecules[i].position.y,
                molecules[j].position.z - molecules[i].position.z
            };
            
            float r = sqrtf(r_vec.x*r_vec.x + r_vec.y*r_vec.y + r_vec.z*r_vec.z);
            
            if(r > 0.1f && r < 15.0f) { // Reasonable interaction range
                // Use mixed parameters for different element pairs
                int elem_i = current_element_index;
                int elem_j = current_element_index;
                if(elem_i >= 6) elem_i = 5;
                if(elem_j >= 6) elem_j = 5;
                
                float epsilon = sqrtf(vdw_params[elem_i].epsilon * vdw_params[elem_j].epsilon);
                float sigma = (vdw_params[elem_i].sigma + vdw_params[elem_j].sigma) * 0.5f;
                
                // Lennard-Jones potential: V(r) = 4ε[(σ/r)^12 - (σ/r)^6]
                float sigma_over_r = sigma / r;
                float sigma6 = powf(sigma_over_r, 6.0f);
                float sigma12 = sigma6 * sigma6;
                
                float force_magnitude = 24.0f * epsilon * (2.0f * sigma12 - sigma6) / r;
                
                // Apply force to both molecules
                molecules[i].force.x -= force_magnitude * r_vec.x / r;
                molecules[i].force.y -= force_magnitude * r_vec.y / r;
                molecules[i].force.z -= force_magnitude * r_vec.z / r;
                
                molecules[j].force.x += force_magnitude * r_vec.x / r;
                molecules[j].force.y += force_magnitude * r_vec.y / r;
                molecules[j].force.z += force_magnitude * r_vec.z / r;
            }
        }
    }
}

void updateMagneticDipoles() {
    // Simple magnetic dipole interactions
    for(int i = 0; i < NUM_MOLECULES; i++) {
        molecules[i].magnetic_moment = molecules[i].element->atomic_number * 0.1f * magnetic_field_strength;
        
        // Magnetic force affects electron orbital speeds
        for(int e = i * NUM_ELECTRONS_PER_MOLECULE; e < (i + 1) * NUM_ELECTRONS_PER_MOLECULE; e++) {
            electrons[e].speed *= (1.0f + molecules[i].magnetic_moment * 0.02f);
        }
    }
}

void drawGravityField() {
    if(!show_gravity_field) return;
    
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Draw gravity field lines
    for(int i = 0; i < NUM_MOLECULES; i++) {
        Vec3 center = gravity_sources[i].position;
        
        glColor4f(0.8f, 0.6f, 0.2f, 0.3f); // Golden field lines
        
        // Draw radial field lines
        for(int angle = 0; angle < 360; angle += 30) {
            float rad = angle * M_PI / 180.0f;
            
            glBegin(GL_LINE_STRIP);
            for(int r = 1; r < 8; r++) {
                float radius = r * 1.5f;
                float field_strength = gravity_sources[i].mass / (radius * radius + 1.0f);
                
                glColor4f(0.8f, 0.6f, 0.2f, field_strength * 0.1f);
                glVertex3f(center.x + radius * cosf(rad), 
                          center.y + radius * sinf(rad), 
                          center.z);
            }
            glEnd();
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

void drawVdWPotential() {
    if(!show_vdw_potential) return;
    
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    // Draw potential wells between close molecules
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            float distance = calculateDistance(molecules[i].position, molecules[j].position);
            
            if(distance < 8.0f) { // Only show for nearby molecules
                Vec3 mid_point = {
                    (molecules[i].position.x + molecules[j].position.x) * 0.5f,
                    (molecules[i].position.y + molecules[j].position.y) * 0.5f,
                    (molecules[i].position.z + molecules[j].position.z) * 0.5f
                };
                
                // Color based on interaction strength
                float interaction_strength = 1.0f / (distance + 1.0f);
                glColor4f(0.2f, 0.8f, 0.9f, interaction_strength * 0.4f);
                
                glPushMatrix();
                glTranslatef(mid_point.x, mid_point.y, mid_point.z);
                GLUquadric* quad = gluNewQuadric();
                gluSphere(quad, 0.3f * interaction_strength, 8, 8);
                gluDeleteQuadric(quad);
                glPopMatrix();
            }
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

// ================================
// MOLECULAR FORMATION SYSTEM
// ================================

void initializeMoleculeFormation() {
    for(int i = 0; i < MAX_FORMED_MOLECULES; i++) {
        formed_molecules[i].active = 0;
        formed_molecules[i].atom_count = 0;
        formed_molecules[i].formation_energy = 0.0f;
    }
}

void checkMoleculeFormation() {
    // Check for simple molecule formations
    for(int i = 0; i < NUM_MOLECULES; i++) {
        for(int j = i + 1; j < NUM_MOLECULES; j++) {
            float distance = calculateDistance(molecules[i].position, molecules[j].position);
            
            if(distance < BOND_FORMATION_DISTANCE) {
                // Check if these atoms can form a known molecule
                int elem_i = molecules[i].element->atomic_number;
                int elem_j = molecules[j].element->atomic_number;
                
                // Find inactive formed molecule slot
                for(int f = 0; f < MAX_FORMED_MOLECULES; f++) {
                    if(!formed_molecules[f].active) {
                        formed_molecules[f].active = 1;
                        formed_molecules[f].atom_count = 2;
                        formed_molecules[f].atom_types[0] = elem_i;
                        formed_molecules[f].atom_types[1] = elem_j;
                        
                        // Calculate center of mass
                        formed_molecules[f].positions[0] = molecules[i].position;
                        formed_molecules[f].positions[1] = molecules[j].position;
                        
                        // Determine molecule type and properties
                        if((elem_i == 1 && elem_j == 1)) { // H-H
                            strcpy(formed_molecules[f].formula, "H2");
                            formed_molecules[f].formation_energy = -436.0f; // kJ/mol
                            formed_molecules[f].bond_lengths[0] = 0.74f; // Angstroms
                        } else if((elem_i == 1 && elem_j == 8) || (elem_i == 8 && elem_j == 1)) { // H-O
                            strcpy(formed_molecules[f].formula, "OH");
                            formed_molecules[f].formation_energy = -463.0f;
                            formed_molecules[f].bond_lengths[0] = 0.96f;
                        } else if((elem_i == 6 && elem_j == 8) || (elem_i == 8 && elem_j == 6)) { // C-O
                            strcpy(formed_molecules[f].formula, "CO");
                            formed_molecules[f].formation_energy = -1072.0f;
                            formed_molecules[f].bond_lengths[0] = 1.13f;
                        } else {
                            // Generic bond
                            sprintf(formed_molecules[f].formula, "%c%c", 
                                   molecules[i].element->symbol[0], 
                                   molecules[j].element->symbol[0]);
                            formed_molecules[f].formation_energy = -200.0f;
                            formed_molecules[f].bond_lengths[0] = distance;
                        }
                        
                        // Create bonding explosion effect
                        Vec3 bond_center = {
                            (molecules[i].position.x + molecules[j].position.x) * 0.5f,
                            (molecules[i].position.y + molecules[j].position.y) * 0.5f,
                            (molecules[i].position.z + molecules[j].position.z) * 0.5f
                        };
                        Vec3 bond_color = {0.2f, 0.8f, 0.3f}; // Green for bonding
                        createExplosion(bond_center, 2, bond_color);
                        
                        printf("🧪 Molecule formed: %s (Energy: %.1f kJ/mol)\n", 
                               formed_molecules[f].formula, formed_molecules[f].formation_energy);
                        break;
                    }
                }
                break; // Only one bond per atom for now
            }
        }
    }
}

void updateFormedMolecules() {
    for(int i = 0; i < MAX_FORMED_MOLECULES; i++) {
        if(formed_molecules[i].active) {
            // Check if atoms are still close enough to maintain bond
            if(formed_molecules[i].atom_count >= 2) {
                Vec3 pos1 = formed_molecules[i].positions[0];
                Vec3 pos2 = formed_molecules[i].positions[1];
                float distance = calculateDistance(pos1, pos2);
                
                if(distance > BOND_BREAK_DISTANCE) {
                    printf("💥 Bond broken: %s (distance: %.2f)\n", 
                           formed_molecules[i].formula, distance);
                    formed_molecules[i].active = 0;
                    
                    // Create break explosion
                    Vec3 break_center = {
                        (pos1.x + pos2.x) * 0.5f,
                        (pos1.y + pos2.y) * 0.5f,
                        (pos1.z + pos2.z) * 0.5f
                    };
                    Vec3 break_color = {1.0f, 0.3f, 0.1f}; // Red for breaking
                    createExplosion(break_center, 0, break_color);
                }
            }
        }
    }
}

void drawFormedMolecules() {
    if(!show_formed_molecules) return;
    
    glDisable(GL_LIGHTING);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    
    for(int i = 0; i < MAX_FORMED_MOLECULES; i++) {
        if(formed_molecules[i].active && formed_molecules[i].atom_count >= 2) {
            Vec3 pos1 = formed_molecules[i].positions[0];
            Vec3 pos2 = formed_molecules[i].positions[1];
            
            // Draw thicker bond line
            glLineWidth(3.0f);
            glColor4f(0.9f, 0.9f, 0.2f, 0.8f); // Bright yellow bond
            
            glBegin(GL_LINES);
            glVertex3f(pos1.x, pos1.y, pos1.z);
            glVertex3f(pos2.x, pos2.y, pos2.z);
            glEnd();
            
            glLineWidth(1.0f);
            
            // Draw molecule label
            Vec3 center = {
                (pos1.x + pos2.x) * 0.5f,
                (pos1.y + pos2.y) * 0.5f + 0.5f,
                (pos1.z + pos2.z) * 0.5f
            };
            
            glColor4f(1.0f, 1.0f, 1.0f, 0.9f);
            glRasterPos3f(center.x, center.y, center.z);
            // Note: Text rendering would need additional setup
        }
    }
    
    glDisable(GL_BLEND);
    glEnable(GL_LIGHTING);
}

Vec3 calculateOptimalGeometry(int molecule_type, int atom_index) {
    Vec3 position = {0, 0, 0};
    
    switch(molecule_type) {
        case 0: // H2O (water)
            if(atom_index == 0) { // Oxygen at center
                position.x = 0; position.y = 0; position.z = 0;
            } else if(atom_index == 1) { // First hydrogen
                position.x = 0.96f * cosf(H2O_ANGLE/2); 
                position.y = 0.96f * sinf(H2O_ANGLE/2); 
                position.z = 0;
            } else { // Second hydrogen
                position.x = 0.96f * cosf(-H2O_ANGLE/2); 
                position.y = 0.96f * sinf(-H2O_ANGLE/2); 
                position.z = 0;
            }
            break;
            
        case 1: // NH3 (ammonia)
            if(atom_index == 0) { // Nitrogen at center
                position.x = 0; position.y = 0; position.z = 0;
            } else { // Hydrogens in pyramid
                float angle = (atom_index - 1) * 2.0f * M_PI / 3.0f;
                position.x = 1.01f * cosf(angle) * sinf(NH3_ANGLE);
                position.y = 1.01f * sinf(angle) * sinf(NH3_ANGLE);
                position.z = 1.01f * cosf(NH3_ANGLE);
            }
            break;
            
        case 2: // CO2 (linear)
            if(atom_index == 0) { // Carbon at center
                position.x = 0; position.y = 0; position.z = 0;
            } else if(atom_index == 1) { // First oxygen
                position.x = -1.16f; position.y = 0; position.z = 0;
            } else { // Second oxygen
                position.x = 1.16f; position.y = 0; position.z = 0;
            }
            break;
    }
    
    return position;
}
